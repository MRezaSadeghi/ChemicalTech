clc; clear; close all;

%% Data
xf = [0.3, 0.7];
xD = [0.95, 0.05];
xB = [0.04, 0.96];

A = [15.9008, 16.0137];
B0 = [2788.51, 3096.52];
C = [-52.36, -53.67];
Tbp = [353.3, 383.8];

P = 760;
q = 1;

cost_fun = @(T) P - sum(antoine_T2P(T, A, B0, C).*xf);
Tf = fzero(cost_fun, 390);
kf = antoine_T2P(Tf, A, B0, C)/P;
af = kf(1)/kf(2);

cost_fun = @(T) P - sum(antoine_T2P(T, A, B0, C).*xB);
TB = fzero(cost_fun, 390);
kB = antoine_T2P(TB, A, B0, C)/P;
aB = kB(1)/kB(2);

cost_fun = @(T) 1/P - sum(xD./antoine_T2P(T, A, B0, C));
TD = fzero(cost_fun, 390);
kD = antoine_T2P(TD, A, B0, C)/P;
aD = kD(1)/kD(2);

a_avg = (af*aD*aB)^(1/3);

xI = xf(1);
yI = a_avg*xI/(1 + (a_avg-1)*xI);
r0 = (xD(1) - yI)/(xD(1) - xI);
Rmin = r0/(1-r0);
R = 1.3*Rmin;


DF = (xf(1) - xB(1))/(xD(1) - xB(1));
BF = 1 - DF;

VpF = DF*R + q + BF;
LpF = VpF - BF;


McCable(a_avg, xD(1), xB(1), xf(1), R, BF, q);
Fenske_corr(a_avg, xD(1), xB(1), R, Rmin);

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