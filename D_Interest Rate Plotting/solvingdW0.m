function f = solvingdW0(x, parameters)
% Creates function handle with transfer rate (tauo1) become the unknown variable x. 
% This function handle requires parameters and mitigation share (m0s).
% The return value is the utility gains for generation 0
% between with and without social contract with unknown variable tauo1.

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

contract   = x(1);
dW0        = CalcCESdW0(contract, parameters);
f          = dW0;
