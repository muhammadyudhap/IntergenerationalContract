function f = rootk1(X, parameters)
% Creates a function handle to solve for k_{t+1}

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

k1       = X(1);

Zemin    = A * exp(-xi * abs(Emin) );
I0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1)*...
    (1 - ( alpha/(1 - alpha) ) * k0^rho * tauo0);
E0       = (1 - delta) * Emin + zeta*k0 - gamma * M0;
Ze0      = A * exp( -xi * abs(E0) );

f        = k1 + k1 * ( beta^(1/theta) * alpha * Ze0 * ( alpha + (1 - alpha)*...
    k1^(-rho) )^(1/rho - 1) * (1 + tauo1) )^( theta / (theta - 1) ) - I0 * (1 - m0);

