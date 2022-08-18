%%
clear;
clc;

% Parameters
alpha    = 1/3;     %parameters(1)
beta     = .7;     %parameters(2)
delta    = 0;       %parameters(3)
gamma    = 1;       %parameters(4)
rho      = 0.0001;   %parameters(5)
theta    = 0;     %parameters(6)
A        = 3;       %parameters(7)
m0       = 0;       %parameters(8) No social contract case
m1       = 0;       %parameters(9)
tauo0    = 0;       %parameters(10)
tauo1    = 0;       %parameters(11) No social contract case
tauo2    = 0;       %parameters(12)
k0       = 2;       %parameters(13)
M0       = 0;       %parameters(14)
Emin     = 0;       %parameters(15)
zeta     = 1;     %parameters(16)
xi       = 1;       %parameters(17)

parameters = [ alpha beta delta gamma rho theta A m0 m1 tauo0 tauo1 tauo2 k0...
    M0 Emin zeta xi];

N1       = 99; 
N2       = 500;
di       = 0.01; 
tauo1s   = 0.01;
T1s      = 0.001;

%% 
guessk1 = 0.5;
 for i = 2:N1
    thetasim1      = (i-1)*di; 
    parameters(6)  = thetasim1;
    parameters(11)  = 0;
    [k1m01(i), fval, eflag, out] = fsolve(@(X) rootk1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
   
    parameters(11)  = tauo1s;
    [k1m11(i), fval, eflag, out] = fsolve(@(X) rootk1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    thetam1(i)     =(i-1)*di; 
end;
diffk11    = k1m11-k1m01;

for i = 2:N2
    thetasim2      = (i-1)*(-di); 
    parameters(6) = thetasim2;
    parameters(11) = 0;
    [k1m02(i), fval, eflag, out] = fsolve(@(X) rootk1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    
    parameters(11)  = tauo1s;
    [k1m12(i), fval, eflag, out] = fsolve(@(X) rootk1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    thetam2(i)     =(i-1)*(-di);
end;
diffk12    = k1m12-k1m02;

diffk1all_sub  = [diffk11, diffk12];


thetaall   = [thetam1, thetam2];
thetaall(thetaall==0)=NaN;
%%
T0       = 0;       %parameters(10)
T1       = 0;       %parameters(11) No social contract case
T2       = 0;       %parameters(12)  

parameters = [ alpha beta delta gamma rho theta A m0 m1 T0 T1 T2 k0...
    M0 Emin zeta xi];

%%
guessk1 = 0.5;
 for i = 2:N1
    thetasim1      = (i-1)*di; 
    parameters(6)  = thetasim1;
    parameters(11)  = 0;
    [k1m01(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));   
    
    parameters(6)  = thetasim1;
    parameters(11)  = T1s;
    [k1m11(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
end;

diffk11_lum    = k1m11-k1m01;

for i = 2:N2
    thetasim2      = (i-1)*(-di); 
    parameters(6) = thetasim2;
    parameters(11) = 0;
    [k1m02(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
  
    parameters(6) = thetasim2;
    parameters(11)  = T1s;
    [k1m12(i), fval, eflag, out] = fsolve(@(X) lumsolvek1(X, parameters),...
        guessk1,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
end;
diffk12_lum    = k1m12-k1m02;

diffk1all_lum  = [diffk11_lum, diffk12_lum];

%% figure
h = figure;
plot(thetaall,diffk1all_sub,'-','LineWidth',2,'Color','blue');
hold on;
plot(thetaall,diffk1all_lum,'--','LineWidth',2,'Color','red');
grid on;
legend({'Capital Subsidy','Lump-sum'},'Interpreter','latex','FontSize',14,'Location','northwest');
xlabel('Intertemporal Consumption Substitution $\theta$','Interpreter','latex','FontSize',14,'Color','black');
ylabel('Changes in Savings $k_{t+1}$','Interpreter','latex','FontSize',14,'Color','black');
%dim = [.4 .7 .2 .2];
%str = {'$\theta=$' theta '$\rho =$' rho '$k0 = $' k0 '$A = $' A '$\tau^o_{t+1}=$' tauo1 '$\mathcal{T}_{t+1}=$' T1 };
%annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');
hold off;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,'theta_k1_dittransfer_both_A15_k1.pdf','-dpdf','-r0');
