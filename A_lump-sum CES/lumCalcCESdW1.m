function f = lumCalcCESdW1(contract, parameters)
% Calculates the the utility gains for generation 1 between 
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
w0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1);
E0       = (1 - delta) * Emin + zeta * k0 - gamma * M0;
Ze0      = A * exp( -xi * abs(E0) );
%% generation t
% solve for k1 under no social contract
guessk1  = 0.5;
k1n      = fsolve(@(X) lumsolvek1(X, parameters), guessk1,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
%k1n denote savings k1 under no social contract

% solve for k1 with existing Social contract
parameters(8) = m0s;
parameters(11)= T1s;
guessk1       = 0.5;
k1s           = fsolve(@(X) lumsolvek1(X, parameters), guessk1,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
% k1s denote savings k1 when social contract exist

%% generation t+1
% No social contract case
% reassigning (M0, T1) = (0,0)
m0            = 0; 
T1            = 0;
parameters(8) = m0;
parameters(11)= T1;
%Solve for k2
parameters(18)= k1n; % add k1 into the parameters for solving k2
guessk2       = 0.5;
k2n           = fsolve(@(X) lumsolvek2(X, parameters), guessk2,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
%k2n denote savings k2 under no social contract

%Calculating other variables
w1       = (1 - alpha) * Ze0 * ( alpha * k1n^rho + (1 - alpha) )^(1/rho - 1);
E1       = (1 - delta) * E0 + zeta * k1n - gamma * m0 * w0;
Ze1      = A * exp( -xi * abs(E1) );
R2       = alpha * Ze1 * ( alpha + (1 - alpha)* k2n^(-rho) )^(1/rho - 1);

cy1      = w1 * (1 - m1) - T1 - k2n;
co2      = R2 * k2n + T2;
U1A      = (cy1^theta)/theta + beta * (co2^theta)/theta; %the utility under no social contract generation t+1

% With social contract case
parameters(8) = m0s;
parameters(11)= T1s; % Trial and error again solving for dW1 = 0
%Solve for k2
parameters(18)  = k1s; % add k1 with social contract into the parameters for solving k2
guessk2         = 0.5;
k2s             = fsolve(@(X) lumsolvek2(X, parameters), guessk2,...
    optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000)); 
%k2s denote savings k2 when social contract exist

%Recalculating other variables
w1      = (1 - alpha) * Ze0 * ( alpha * k1s^rho + (1 - alpha) )^(1/rho - 1);
E1      = (1 - delta) * E0 + zeta*k1s - gamma * m0s * w0;
Ze1     = A * exp( -xi * abs(E1) );
R2      = alpha * Ze1 * ( alpha + (1 - alpha)* k2s^(-rho) )^(1/rho - 1);

cy1     = w1 * (1 - m1) - T1s - k2s;
co2     = R2 * k2s + T2;
U1B     = (cy1^theta)/theta + beta * (co2^theta)/theta; %the utility with social contract generation t+1
dW1     = U1B - U1A;  %The goal is to search for T1 when dW1 near to zero
f       = dW1;
