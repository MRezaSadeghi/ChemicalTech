clc; clear; close all;

F = 1;
alpha = 0.25;
zf = [0.071, 0.215, 0.357, 0.286, 0.071];
A = [15.726, 15.5381, 15.6782, 15.6338, 15.8333];
B = [1872.45, 2032.73, 2154.9, 2348.67, 2477.07];
C = [-25.16, -33.15, -34.42, -40.05, -39.94];
Tc = [369.9, 408.1, 425.2, 460.35, 469.6];

P1 = 765.0;
P1 = P1/0.133322;

Tguess = 340;
cost_fun = @(T1) Rachford_t(T1, P1, A, B, C, alpha, zf);
T1 = fzero(cost_fun, Tguess);
k1 = antoine_T2P(T1, A, B, C)/P1;

x1 = zf./(1 + alpha*(k1 - 1));
y1 = x1.*k1;

V1 = alpha*F;
L1 = F - V1;

T2 = 325;
P2 = 400;
P2 = P2/0.133322;

Pev2 = antoine_T2P(T2, A, B, C);
k2  = Pev2/P2;

alpha_guess = 0.5;
cost_fun = @(a) Rachford_a(k2, a, x1);
alpha2 = fzero(cost_fun, alpha_guess);

x2 = x1./(1 + alpha2.*(k2 - 1));
y2 = k2.*x2;

V2 = alpha2*L1;
L2 = L1 - V2;


Hev = [18.80, 21.30, 22.40, 25.43, 25.80]; %[kJ/mol]
Tbp = [231.1, 261.3, 272.7, 302.6, 309.2];
cp = [63.05, 88.39, 92.07, 122, 129]/1000; %[kJ/mol]


Hev1 = HevCorrection(Hev, Tbp, T1, Tc);
Hev2 = HevCorrection(Hev, Tbp, T2, Tc);

H1L = L1 * sum(x1.*(-Hev1 + cp.*(T1 - T2)));
H2L = L2 * sum(x2.*(-Hev2));


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

function Nmin = Fenske_corr(a, xD, xB)
    Nmin = log((xD/(1-xD))*(1-xB)/xB);
    Nmin = Nmin/log(a);
end

function Ns = McCable(alpha, xD, xF, R)
    a = 1/(alpha-1);
    b = xD/R + alpha*(1+R)/(1-alpha)/R;
    c = xD/R*a;

    delta_1 = (-(a+b) + sqrt((a+b)^2-4*c))/2;
    delta_2 = (-(a+b) - sqrt((a+b)^2-4*c))/2;

    delta = delta_1;
    A = 1/(xD - delta) + 1/(2*delta+a+b);
    B = 1/(xF - delta) + 1/(2*delta+a+b);
    Ns = log(A/B)/log(-(a+delta)/(b+delta));
    
end