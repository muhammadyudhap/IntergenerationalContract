function f = lumsolvedW0(x, parameters, m0s)
% Creates function handle with lump sum transfer T1s become the unknown variable x. 
% This function handle requires parameters and mitigation effort M0s.
% The return value is the utility gains for generation 0
% between with and without social contract with unknown variable T1s.

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

contract   = [m0s x(1)];
dW0        = CalcLogdW0(contract, parameters);
f          = dW0;
