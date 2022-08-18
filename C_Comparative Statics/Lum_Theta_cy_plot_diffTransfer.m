%%
clear;
clc;

% Parameters
alpha    = 1/3;     %parameters(1)
beta     = .7;      %parameters(2)
delta    = 0;       %parameters(3)
gamma    = 1;       %parameters(4)
rho      = 0.0001;      %parameters(5)
theta    = 0;       %parameters(6)
A        = 3;       %parameters(7)
m0       = 0;       %parameters(8) No social contract case
m1       = 0;       %parameters(9)
T0       = 0;       %parameters(10)
T1       = 0;       %parameters(11) No social contract case
T2       = 0;       %parameters(12)
k0       = 2;       %parameters(13)
M0       = 0;       %parameters(14)
Emin     = 0;       %parameters(15)
zeta     = 1;       %parameters(16)
xi       = 1;       %parameters(17)

parameters = [ alpha beta delta gamma rho theta A m0 m1 T0 T1 T2 k0...
    M0 Emin zeta xi];

N1       = 100; 
N2       = 500;
di       = 0.01; 
Mos      = 0.01;
T1s      = 0.001;

Zemin    = A * exp( -xi * abs(Emin) );
w0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1);
E0       = (1 - delta) * Emin + zeta * k0 - gamma * M0;
Ze0      = A * exp( -xi * abs(E0) );

%%
guessk1 = 0.5;
 for i = 2:N1
    thetasim1      = (i-1)*di; 
    parameters(6)  = thetasim1;
    parameters(11)  = 0;
    [k1m01(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));   
    R100(i)      = alpha * Ze0 * ( alpha + (1 - alpha)* k1m01(i)^(-rho) )^(1/rho - 1);
    cy00(i)      = w0 * (1 - parameters(8)) - T0 - k1m01(i);
    co00(i)      = R100(i) * k1m01(i) + parameters(11);
    
    parameters(6)  = thetasim1;
    parameters(11)  = T1s;
    [k1m11(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    R101(i)      = alpha * Ze0 * ( alpha + (1 - alpha)* k1m11(i)^(-rho) )^(1/rho - 1);
    cy01(i)      = w0 * (1 - parameters(8)) - T0 - k1m11(i);
    co01(i)      = R101(i) * k1m11(i) + parameters(11);
    thetam1(i)    =(i-1)*di; 
end;
diffk11    = k1m11-k1m01;
diffcy01   = cy01-cy00;
diffco11   = co01-co00;
%diffk11(diffk11==0)= NaN;

%%
for i = 2:N2
    thetasim2      = (i-1)*(-di); 
    parameters(6) = thetasim2;
    parameters(11) = 0;
    [k1m02(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    R110(i)      = alpha * Ze0 * ( alpha + (1 - alpha)* k1m02(i)^(-rho) )^(1/rho - 1);
    cy10(i)      = w0 * (1 - parameters(8)) - T0 - k1m02(i);
    co10(i)      = R110(i) * k1m02(i) + parameters(11);
    
    
    parameters(6) = thetasim2;
    parameters(11)  = T1s;
    [k1m12(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    R111(i)      = alpha * Ze0 * ( alpha + (1 - alpha)* k1m12(i)^(-rho) )^(1/rho - 1);
    cy11(i)      = w0 * (1 - parameters(8)) - T0 - k1m12(i);
    co11(i)      = R111(i) * k1m12(i) + parameters(11);
    thetam2(i)   =(i-1)*(-di);
end;
diffk12    = k1m12-k1m02;
diffcy02   = cy11-cy10;
diffco12   = co11-co10;
%diffk12(diffk12==0)= NaN;


%%
diffk1all  = [diffk11, diffk12];
thetaall   = [thetam1, thetam2];
diffcyall  = [diffcy01, diffcy02];
diffcoall  = [diffco11,diffco12];
thetaall(thetaall==0)=NaN;


h = figure;
plot(thetaall,diffcyall,'-','LineWidth',2,'Color','blue');
hold on;
plot(thetaall,diffcoall,'--','LineWidth',2,'Color','red');
plot(thetaall,diffk1all,':','LineWidth',2,'Color','green');
grid on;
legend({'$\Delta c^y_t$','$\Delta c^o_{t+1}$', '$\Delta k_{t+1}$'},'Interpreter','latex','FontSize',18,'Location','northwest');
xlabel('Intertemporal Substitution Parameter\hspace{0.15cm} $\theta$','Interpreter','latex','FontSize',14,'Color','black');
ylabel('Changes in Consumption/Savings','Interpreter','latex','FontSize',14,'Color','black');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,'CompStatLum_T1_rhomin3.pdf','-dpdf','-r0');
