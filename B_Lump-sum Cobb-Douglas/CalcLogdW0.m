function f = CalcLogdW0(contract, parameters)
% Calculates the the utility gains for generation 0 between 
% with and without lump sum social contract.
% The lump sum social contract (M0, T1) must be assigned first.

% Parameters
alpha      = parameters(1);
beta       = parameters(2);
delta      = parameters(3);
gamma      = parameters(4);
rho        = parameters(5);
theta      = parameters(6);
A          = parameters(7);
m0         = parameters(8);
m1         = parameters(9);
T0         = parameters(10);
T1         = parameters(11);
T2         = parameters(12);
k0         = parameters(13);
M0         = parameters(14);
Emin       = parameters(15);
zeta       = parameters(16);
xi         = parameters(17);

parameters = [ alpha beta delta gamma rho theta A m0 m1 T0 T1 T2 k0...
    M0 Emin zeta xi];
%% The contract
m0s       = contract(1);
T1s       = contract(2); 


Zemin    = A * exp( -xi * abs(Emin) );
w0       = (1 - alpha) * Zemin * k0^alpha;
E0       = (1 - delta) * Emin + zeta * k0 - gamma * M0;
Ze0      = A * exp( -xi * abs(E0) );

%% generation t
% solve for k1 under no social contract
guessk1  = 0.5;
k1n      = fsolve(@(X) logsolvek1(X, parameters), guessk1,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
%k1n denote savings k1 under no social contract

%calculating other variables under no social contract
R1       = alpha * Ze0 * k1n^(alpha - 1);
cy0      = w0 * (1 - m0) - T0 - k1n;
co1      = R1 * k1n + T1;
U0A      = log(cy0) + beta * log(co1); %the utility under no social contract

%% Social contract

%solve for k1, with assigned contract (M0S, tauo1)
parameters(8) = m0s;
parameters(11)= T1s;
guessk1       = 0.5;
k1s           = fsolve(@(X) logsolvek1(X, parameters), guessk1,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
%k1s denote savings k1 when social contract exist

%Recalculating other variables with social contract
R1      = alpha * Ze0 * k1s^(alpha - 1);
cy0     = w0 * (1 - m0s) - T0 - k1s;
co1     = R1 * k1s + T1s;

U0B     = log(cy0) + beta * log(co1); %the utility under existing social contract
dW0     = U0B - U0A; %The goal is to search for T1 when dW0 near to zero
f       = dW0;
