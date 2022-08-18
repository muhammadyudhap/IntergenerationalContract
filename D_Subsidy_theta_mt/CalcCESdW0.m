function f = CalcCESdW0(contract, parameters)
% Calculates the the utility gains for generation 0 between 
% with and without social contract.
% The social contract (m0, tauo1) must be assigned first.

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
tauo0      = parameters(10);
tauo1      = parameters(11);
tauo2      = parameters(12);
k0         = parameters(13);
M0         = parameters(14);
Emin       = parameters(15);
zeta       = parameters(16);
xi         = parameters(17);

parameters = [ alpha beta delta gamma rho theta A m0 m1 tauo0 tauo1 tauo2 k0...
    M0 Emin zeta xi];
%% The contract
theta  = contract(1);
m0s    = contract(2); 


Zemin    = A * exp( -xi * abs(Emin) );
I0       = (1 - alpha) * Zemin * k0^alpha *(1 - ( alpha/(1 - alpha) )* tauo0);
E0       = (1 - delta) * Emin + zeta * k0 - gamma * M0;
Ze0      = A * exp( -xi * abs(E0) );

%% generation t
% solve for k1 under no social contract
parameters(11) = 0; 
guessk1  = 0.5;
k1n      = fsolve(@(X) rootk1(X, parameters), guessk1,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
%k1n denote savings k1 under no social contract

%calculating other variables under no social contract
R1       = alpha * Ze0 * k1n^(alpha-1);
cy0      = I0 * (1 - m0)  - k1n;
co1      = R1 * (1 + 0) * k1n;
U0A      = (cy0^theta)/theta + beta * (co1^theta)/theta; %the utility under no social contract

%% Social contract

%solve for k1, with assigned contract (M0S, tauo1)
parameters(8) = m0s;
parameters(11)= tauo1;
guessk1       = 0.5;
k1s           = fsolve(@(X) rootk1(X, parameters), guessk1,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
%k1s denote savings k1 when social contract exist

%Recalculating other variables with social contract
R1      = alpha * Ze0 * k1s^(alpha-1);
cy0     = I0 * (1 - m0s) - k1s;
co1     = R1 * (1 + tauo1) * k1s;
U0B     = (cy0^theta)/theta + beta * (co1^theta)/theta; %the utility under existing social contract
dW0     = U0B - U0A; %The goal is to search for tauo1 when dW0 near to zero
f       = dW0;
