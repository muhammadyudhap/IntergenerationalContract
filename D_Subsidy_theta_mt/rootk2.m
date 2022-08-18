function f = rootk2(X, parameters)
% Creates a function handle to solve for k_{t+2}

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
tauo0    = parameters(10);
tauo1    = parameters(11);
tauo2    = parameters(12);

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
I0       = (1 - alpha) * Zemin * k0^alpha *(1 - ( alpha/(1 - alpha) )* tauo0);
I1       = (1 - alpha) * Ze0 * k1^alpha *(1 - ( alpha/(1 - alpha) )* tauo1);
E1       = (1 - delta) * E0 + zeta*k1 - gamma * m0 * I0;
Ze1      = A * exp( -xi * abs(E1) );

f = k2 + k2 * ( beta^(1/theta) * alpha * Ze1 * k2^(alpha-1)*...
    (1 + tauo2))^( theta / (theta - 1) ) - I1 * (1 - m1);

