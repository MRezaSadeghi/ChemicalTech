clc; clear; close all;

% inline finctoins
a2A = @(a) a/(1-a);
A2a = @(A) A/(1+A);


y0 = 5e-2;
yN = 5e-3;
Q = 1500; %m3/h
G0 = 1500/22.414;
G = G0*(1 - y0);
P0 = 1;

Y0 = a2A(y0);
YN = a2A(yN);

w = 0.2/100;

mu_o = 240.475;
mu_b = 78.114; 
mu_mix = 1/(w/mu_o + (1-w)/mu_b);
xN1 = w/mu_o * mu_mix;

A = [15.9008, 16.1510];
B = [2788.51, 4294.55];
C = -[52.36, 124];

T = 300;

P = antoine_T2P(T, A, B, C);
Pev_b = P(1)/760;
Pev_o = P(2)/760;

XN1 = a2A(xN1);

K = Pev_b/P0;

[LG, ytg] = riccati_abs_solver(K, yN, xN1, Y0);
LG = 1.5*LG;

i = 0;
value = Y0;
while value >= yN
    i = i+1;
    value = riccati_abs(G, L, K, yN, xN1, value);
    fprintf("i = %d: Yi=%2.5f\n", [i, value]);
    
end


%% Funcitons

function [LG, ytg]= riccati_abs_solver(K, yN, xN1, Y0)
    syms LG;
    
    B = K/(K-1);
    A = LG/(K-1) + (yN - LG*xN1);
    D = B*(yN - LG*xN1);
    delta = (A - B)^2  + 4*D;

    S = solve(delta, LG);
    S = double(S);

    LG = S(1);
    ytg1 = double((subs(A, LG)- subs(B, LG))/2);
    fprintf("L/G = %2.4f and Ytg = %2.4f\n", [LG, ytg1])

    LG = S(2);
    ytg2 = double((subs(A, LG)- subs(B, LG))/2);
    fprintf("L/G = %2.4f and Ytg = %2.4f\n", [LG, ytg2])

    ytg = [ytg1, ytg2];
    ytg = ytg(ytg > yN);
    ytg = ytg(ytg < Y0);
    LG = S([ytg1, ytg2]==ytg);

end

% Antoine conver P to T
function P = antoine_T2P(T, A, B, C)
    lnP = A - B./(T+C);
    P = exp(lnP);
end

function N = Kremser_findN_ybased(E, y0, yN, yeq)
    A = yeq - yN;
    B = (yeq - y0) - E*(yN - y0);
    N = log(A/B)/log(E);
    fprintf("N of Stages: %2.2f = %d\n", [N, ceil(N)])
    N = ceil(N);
end

function yn = Kremser_find_yn_ybased(E, n, N, y0, yeq)
    yn = -(E^(N+1) - E^n)/(E^(N+1) - 1) * (yeq - y0);
    yn = yn + yeq;
end

function N = Kremser_findN_xbased(E, y0, yN, yeq)
    E = 1/E;
    A = yeq - yN;
    B = (yeq - y0) - E*(yN - y0);
    N = log(A/B)/log(E);
end

function yn = Kremser_find_xn_ybased(E, n, y0, yeq)
    E = 1/E;
    yn = -(E^(N+1) - E^n)/(E^(N+1) - 1) * (yeq - y0);
    yn = yn + yeq;
end

function yn1 = riccati_abs(G, L, K, yN, xN1, yn)
    B = K/(K-1);
    A = L/G/(K-1) + (yN - L*xN1/G);
    D = B*(yN - L*xN1/G);
    yn1 = (D - B*yn)/(yn - A);

end

function Yn = Xn2Yn(Xn, K)
    Yn = K.*Xn./(1 + (1-K).*Xn);
end

function Xn = Yn2Xn(Yn, K)
    Xn = Yn./(K - (1-K).*Yn);
end
