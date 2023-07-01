clc; clear; close all;

% inline finctoins
a2A = @(a) a/(1-a);
A2a = @(A) A/(1+A);


% data
y0 = 0.003;
xN1 = 0;
yN = 0.0001;

YN = a2A(yN);
Y0 = a2A(y0);
XN1 = a2A(xN1);

H = 500;
P = 1;
k = H/P;

x1 = y0/k;
X1 = a2A(x1);

LG = (YN - Y0)/(XN1 - X1);
LG = 1.2*LG;

E = k/LG;

yeq = k*XN1;

N = Kremser_findN_ybased(E, Y0, YN, yeq)

Kremser_find_yn_ybased(E, N, N, y0, yeq)

%% Funcitons

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

function [A, B, D, yn1] = riccati_abs(G, L, K, yN, xN1, yn)
    B = K/(K-1);
    A = L/G/(K-1) + (yN - L*xN1/G);
    D = B*(yN - L*xN1/G);
    yn1 = (D - B*yn)/(yn - A);

    if K<1
        fprintf("K in lower than ONE, tangency check\n")
        delta = (A - B)^2  - 4*D;
        if delta >= 0
            fprinft("There is a coincident with delta %2.2f", delta)
        end
    end
end

function Yn = Xn2Yn(Xn, K)
    Yn = K.*Xn./(1 + (1-K).*Xn);
end

function Xn = Yn2Xn(Yn, K)
    Xn = Yn./(K - (1-K).*Yn);
end
