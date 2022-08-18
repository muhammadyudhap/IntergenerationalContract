function f = lumsolvek2(X, parameters)

% Parameters
alpha    = parameters(1);
beta     = parameters(2);
delta    = parameters(3);
gamma    = parameters(4);
rho      = parameters(5);
theta    = parameters(6);
A        = parameters(7);

m0       = parameters(8);
m1       = parameters(9);
T0       = parameters(10);
T1       = parameters(11);
T2       = parameters(12);

k0       = parameters(13);
M0       = parameters(14);
Emin     = parameters(15);
zeta     = parameters(16);
xi       = parameters(17);

k1       = parameters(18);

k2       = X(1);

E0       = (1 - delta) * Emin + zeta*k0 - gamma * M0;
Zemin    = A * exp( -xi * abs(Emin) );
Ze0      = A * exp( -xi * abs(E0) );
w0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1);
w1       = (1 - alpha) * Ze0 * ( alpha * k1^rho + (1 - alpha) )^(1/rho - 1);
E1       = (1 - delta) * E0 + zeta*k1 - gamma * m0 * w0;
Ze1      = A * exp( -xi * abs(E1) );

f        = k2 + k2 * ( beta^(1/theta) * alpha * Ze1 * ( alpha + (1 - alpha) * ...
    k2^(-rho) )^(1/rho - 1) )^( theta / (theta - 1) ) - w1 * (1-m1) + T1 + T2 * ...
    (beta * alpha * Ze1 * ( alpha + (1 - alpha) * k2^(-rho) )^(1/rho - 1) )^( 1 / (theta - 1) );

