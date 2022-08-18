function [U0, U1, k1, k2] = welfarelog(parameters);

% Parameters
alpha    = parameters(1);
beta     = parameters(2);
delta    = parameters(3);
gamma    = parameters(4);
xi       = parameters(5);

m0       = parameters(6);
m1       = parameters(7);
tauo0    = parameters(8);
tauo1    = parameters(9);
tauo2    = parameters(10);

k0       = parameters(11);
M0       = parameters(12);
Emin     = parameters(13);
zeta     = parameters(14);
A        = parameters(15);

I0  = A * exp(-xi * Emin) * ( 1 - alpha ) * k0^alpha * ( 1 - alpha / ( 1 - alpha ) * tauo0 );
k1  = beta / ( 1 + beta ) * ( 1 - m0) * I0;
E0  = ( 1 - delta ) * Emin + zeta*k0 - gamma * M0;
cy0 = 1 / ( 1 + beta ) * ( 1 - m0 ) * I0;
co1 = beta / ( 1 - beta ) * A * exp(-xi * E0) * alpha * k1 ^ ( alpha - 1 ) * ( 1 + tauo1 ) * ( 1 - m0 ) * I0;

I1  = A * exp(-xi * E0)   * ( 1 - alpha ) * k1^alpha * ( 1 - alpha / ( 1 - alpha ) * tauo1 );
k2  = beta / ( 1 + beta ) * ( 1 - m1) * I1;
E1  = ( 1 - delta ) * E0 + zeta*k1 - gamma * m0 * I0;
cy1 = 1 / ( 1 + beta ) * ( 1 - m1 ) * I1;
co2 = beta / ( 1 - beta ) * A * exp(-xi * E1) * alpha * k2 ^ ( alpha - 1 ) * ( 1 + tauo2 ) * ( 1 - m1 ) * I1;

U0 = log( cy0 ) + beta * log ( co1 );
U1 = log( cy1 ) + beta * log ( co2 );

