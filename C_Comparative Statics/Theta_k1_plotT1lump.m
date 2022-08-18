%%
clear;
clc;

% Parameters
alpha    = 1/3;     %parameters(1)
beta     = .7;     %parameters(2)
delta    = 0;       %parameters(3)
gamma    = 1;       %parameters(4)
rho      = -2;     %parameters(5)
theta    = 0;     %parameters(6)
A        = 10;       %parameters(7)
m0       = 0;       %parameters(8) No social contract case
m1       = 0;       %parameters(9)
T0       = 0;       %parameters(10)
T1       = 0;       %parameters(11) No social contract case
T2       = 0;       %parameters(12)
k0       = 0.1;       %parameters(13)
M0       = 0;       %parameters(14)
Emin     = 0;       %parameters(15)
zeta     = 1;     %parameters(16)
xi       = 1;       %parameters(17)

parameters = [ alpha beta delta gamma rho theta A m0 m1 T0 T1 T2 k0...
    M0 Emin zeta xi];

N1       = 100; 
N2       = 500;
di       = 0.01; 
T1s      = 0.1;
k1 = fsolve(@(X) lumsolvek1(X, parameters), 5);
Zemin    = A * exp(-xi * abs(Emin) );
w0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1);
E0       = (1 - delta) * Emin + zeta*k0 - gamma * M0;
Ze0      = A * exp( -xi * abs(E0) );
R1       = alpha * Ze0 * ( alpha + (1 - alpha)* k1^(-rho) )^(1/rho - 1);
%%
guessk1 = 0.5;
 for i = 2:N1
    thetasim1      = (i-1)*di; 
    parameters(6)  = thetasim1;
    parameters(11)  = 0;
    k1m01(i) = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    parameters(11)  = T1s;
    k1m11(i) = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    thetam1(i)     =(i-1)*di; 
end;
diffk11    = k1m11-k1m01;
diffk11(diffk11==0)= NaN;
%%
for i = 2:N2
    thetasim2      = (i-1)*(-di); 
    parameters(6) = thetasim2;
    parameters(11) = 0;
    k1m02(i) = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    parameters(11)  = T1s;
    k1m12(i) = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    thetam2(i)     =(i-1)*(-di);
end;
diffk12    = k1m12-k1m02;
diffk12(diffk12==0)= NaN;
%%
diffk1all  = [diffk11, diffk12];
thetaall   = [thetam1, thetam2];
h = figure;
plot(thetaall,diffk1all,'-','LineWidth',2,'Color','blue');
grid on;
xlabel('Intertemporal Substitution \hspace{0.15cm} $\theta$','Interpreter','latex','FontSize',14,'Color','black');
ylabel('Changes in Savings \hspace{0.15cm} $\Delta k_{t+1}$','Interpreter','latex','FontSize',14,'Color','black');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,'SubCES_k2_A2_thetamin01_rhomin3.pdf','-dpdf','-r0');
