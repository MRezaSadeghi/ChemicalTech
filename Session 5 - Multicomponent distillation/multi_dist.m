clc; clear; close all;

%% Data
xf = [0.04, 0.4, 0.5, 0.06];
Aa = [15.8333, 15.8366, 15.8737, 15.9426];
Ba = [2477.07, 2697.55, 2911.32, 3120.29];
Ca = -[39.94, 48.78, 56.51, 63.63];

Tbp = [309.2, 341.9, 371.6, 398.8];
Hev = [6160, 6896, 7576, 8225];
cpv = [31.02, 38.86, 45.15, 52.23];
cpl = [37.614, 51.700, 56.12, 63.89];
Tc = [469.6, 507.4, 540.2, 568.8];

n = length(xf);
P = 760;
q = 1;
F = 100;

syms xD [1, n]
syms xB [1, n]
syms B D


eqns = [F*xf(1) == D*xD(1),...
        F*xf(2) == D*xD(2) + B*xB(2),...
        F*xf(3) == D*xD(3) + B*xB(3),...
        F*xf(4) == B*xB(4),...
        D*xD(2) == 0.98*F*xf(2),...
        D*xD(3) == 0.01*F*xf(3),...
        sum(xD) - 1,...
        sum(xB) - 1,...
        xB(1),...
        xD(4)];

S = solve(eqns);
xD = [S.xD1, S.xD2, S.xD3, S.xD4];
xD = double(xD);
xB = [S.xB1, S.xB2, S.xB3, S.xB4];
xB = double(xB);
B = double(S.B);
D = double(S.D);

lk = 2;
hk = 3;

cost_fun = @(T) P - sum(antoine_T2P(T, Aa, Ba, Ca).*xf);
Tf = fzero(cost_fun, 390);
kf = antoine_T2P(Tf, Aa, Ba, Ca)/P;
af = kf/kf(hk);

cost_fun = @(T) P - sum(antoine_T2P(T, Aa, Ba, Ca).*xB);
TB = fzero(cost_fun, 390);
kB = antoine_T2P(TB, Aa, Ba, Ca)/P;
aB = kB/kB(hk);

cost_fun = @(T) 1/P - sum(xD./antoine_T2P(T, Aa, Ba, Ca));
TD = fzero(cost_fun, 390);
kD = antoine_T2P(TD, Aa, Ba, Ca)/P;
aD = kD/kD(hk);

cost_fun = @(T) P - sum(antoine_T2P(T, Aa, Ba, Ca).*xD);
TDb = fzero(cost_fun, 390);


a_avg = (af.*aD.*aB).^(1/3);
a_lk = a_avg(lk);

Nmin = Fenske_corr_multi(a_lk, xD, xB, lk, hk);

syms theta
S = solve(sum(a_avg.*xf./(a_avg - theta)) - 1 + q);
S = double(S);
theta = S(logical((S>1) .* (S<a_lk)));

Rmin = sum(a_avg.*xD./(a_avg - theta)) - 1;
R = 1.5*Rmin;

phi_F_calc(Nmin, R, Rmin)

V = D*(1+R);
L = D*R;

% HLF = 0;
% HLD = sum(xD.*cpl.*(TDb - Tf));
% HLB = sum(xB.*cpl.*(TB - Tf));
% HVD = sum(xD.*(cpl.*(Tbp - Tf) + Hev + cpv.*(TD - Tbp)));

HVD = sum(xD.*(cpl.*(Tbp - TDb) + Hev + cpv.*(TD - Tbp)));
Qc = V*HVD


HLF = 0;
HLD = sum(xD.*cpl.*(TDb - Tf));
HLB = sum(xB.*cpl.*(TB - Tf));
Qr = Qc - F*HLF + D*HLD + B*HLB
Qr*4.184




%% Functions

% Hev Correction (0.38)
function H2 = HevCorrection(H1, T1, T2, Tc)
    Tc1 = T1./Tc;
    Tc2 = T2./Tc;
    H2 = H1.*((1-Tc2)./(1-Tc1)).^0.38;

end
% Rachford-Rice (with varibale T)
function obj = Rachford_t(T, P, A, B, C, alpha, zf)
    k = antoine_T2P(T, A, B, C)/P;
    obj = sum(zf.*(k-1)./(1 + alpha.*(k-1)));    
end

% Rachford-Rice (with varibale alpha)
function obj = Rachford_a(k, alpha, zf)
    obj = sum(zf.*(k-1)./(1 + alpha.*(k-1)));    
end

% Antoine conver T to P
function T = antoine_P2T(P, A, B, C)
    T = B./(A - log(P)) - C;
end

% Antoine conver P to T
function P = antoine_T2P(T, A, B, C)
    lnP = A - B./(T+C);
    P = exp(lnP);
end

function Nmin = Fenske_corr_multi(a, xD, xB, lk, hk)
    fprintf("Fenske Short Cut (FUG) Method\n")
    Nmin = log((xD(lk)/xB(lk))*(xB(hk)/xD(hk)));
    
    Nmin = Nmin/log(a);
    Nmin = ceil(Nmin);
    fprintf("Nmin stage = %d\n", Nmin)
end

function phi_F_calc(Nmin, R, Rmin)

    FR = (R - Rmin)/(1 + R);
    phi_ed = 0.75 - 0.75*FR^0.5668;
    phi_mo = 1 - exp((1 + 54.4*FR)*(FR - 1)/(11+117.2*FR)/(sqrt(FR)));

    N_ed = ceil((phi_ed + Nmin)/(1 - phi_ed));
    N_mo = ceil((phi_mo + Nmin)/(1 - phi_mo));

    fprintf("Nstage Molokanov = %d\n", N_mo);
    fprintf("Nstage Eduljee = %d\n", N_ed);
    fprintf("Ntray Molokanov = %d\n", N_mo-1);
    fprintf("Ntray Eduljee = %d\n\n", N_ed-1);

end

function Fenske_corr(a, xD, xB, R, Rmin)
    fprintf("Fenske-Gilliland Method\n")
    Nmin = log((xD/(1-xD))*(1-xB)/xB);
    Nmin = Nmin/log(a);
    FR = (R - Rmin)/(1 + R);
    phi_ed = 0.75 - 0.75*FR^0.5668;
    phi_mo = 1 - exp((1 + 54.4*FR)*(FR - 1)/(11+117.2*FR)/(sqrt(FR)));

    N_ed = ceil((phi_ed + Nmin)/(1 - phi_ed));
    N_mo = ceil((phi_mo + Nmin)/(1 - phi_mo));

    fprintf("Nstage Molokanov = %d\n", N_mo);
    fprintf("Nstage Eduljee = %d\n\n", N_ed);
end

function McCable(alpha, xD, xB, xF, R, BF, q)
    fprintf("McCable Method\n")
    a = 1/(alpha-1);
    b = xD/R + alpha*(1+R)/(1-alpha)/R;
    c = xD/R*a;

    delta_1 = (-(a+b) + sqrt((a+b)^2-4*c))/2;
    delta_2 = (-(a+b) - sqrt((a+b)^2-4*c))/2;

    delta = delta_1;
    Avar = 1/(xD - delta) + 1/(2*delta+a+b);
    Bvar = 1/(xF - delta) + 1/(2*delta+a+b);
    Cvar = -(a+delta)/(b+delta);
    Ns = log(Avar/Bvar)/log(Cvar);
    Ns = ceil(Ns);
    fprintf("Ns = %d\n", Ns)

    xfm = (Avar/Cvar^Ns - 1/(2*delta+a+b))^-1 + delta;

    DF = 1-BF;
    bp = -(BF*xB + (R*DF + q - BF)*alpha/(alpha-1))/(R*DF + q);
    cp = -BF*xB/(R*DF + q)/(alpha-1);

    deltap_1 = (-(a+bp) + sqrt((a+bp)^2-4*cp))/2;
    deltap_2 = (-(a+bp) - sqrt((a+bp)^2-4*cp))/2;

    deltap = deltap_1;
    Avar = 1/(xB - deltap) + 1/(2*deltap+a+bp);
    Bvar = 1/(xfm - deltap) + 1/(2*deltap+a+bp);
    Cvar = -(a+deltap)/(bp+deltap);
    Ninf = log(Bvar/Avar)/log(Cvar);
    Ninf = ceil(Ninf);
    fprintf("Ninf = %d\n", Ninf)

    N = Ns + Ninf;
    fprintf("Nstage = %d\n", N);
    fprintf("Ntray = %d\n\n", N-1);
end