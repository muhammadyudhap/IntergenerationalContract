%%
clear;
clc;

% Parameters
alpha    = 1/3;     %parameters(1)
beta     = .7;     %parameters(2)
delta    = 0;       %parameters(3)
gamma    = 1;       %parameters(4)
rho      = 0.0001;     %parameters(5)
theta    = 0.001;     %parameters(6)
A        = 3;       %parameters(7)
m0       = 0;       %parameters(8) No social contract case
m1       = 0;       %parameters(9)
tauo0    = 0;       %parameters(10)
tauo1    = 0;       %parameters(11) No social contract case
tauo2    = 0;       %parameters(12)
k0       = 0;       %parameters(13)
M0       = 0;       %parameters(14)
Emin     = 0;       %parameters(15)
zeta     = 1;     %parameters(16)
xi       = 1;       %parameters(17)

parameters = [ alpha beta delta gamma rho theta A m0 m1 tauo0 tauo1 tauo2 k0...
    M0 Emin zeta xi];

N        = 50; % Number of iterations 
di       = 0.05; % stepsize 
tauo1s   = 0.01;
T1s      = 0.005;

%% 
 guesstmt  = 0.01;
for i = 1:N
    parameters(13) = i*di; %simulating k0
    parameters(11) = tauo1s; 
    [mtvec(i), fval, eflag, out] = fsolve(@(x) solvingdW0(x, parameters),...
        guesstmt,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    guesstmt = mtvec(i);
    parameters(11) = 0; 
    k0       = parameters(13);
    k1       = fsolve(@(X) rootk1(X, parameters), 5);
    Zemin    = A * exp( -xi * abs(Emin) );
    I0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1)*(1 - ( alpha/(1 - alpha) ) * k0^rho * tauo0);
    E0       = (1 - delta) * Emin + zeta * k0 - gamma * M0;
    Ze0      = A * exp( -xi * abs(E0) );
    R1vec(i) = alpha * Ze0 * ( alpha + (1 - alpha)* k1^(-rho) )^(1/rho - 1);
end;

%%
% Parameters
T0       = 0;       %parameters(10) 
T1       = 0;       %parameters(11) 
T2       = 0;       %parameters(12)
parameters = [ alpha beta delta gamma rho theta A m0 m1 T0 T1 T2 k0...
    M0 Emin zeta xi];
for i = 1:N
    parameters(13) = i*di;
    parameters(11) = T1s; 
    [mtvec2(i), fval, eflag, out] = fsolve(@(x) lumsolvedW0(x, parameters),...
        guesstmt,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    guesstmt = mtvec2(i);
    k0       = parameters(13);
    parameters(11) = 0; 
    k1       = fsolve(@(X) lumsolvek1(X, parameters), 5);
    Zemin    = A * exp(-xi * abs(Emin) );
    w0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1);
    E0       = (1 - delta) * Emin + zeta*k0 - gamma * M0;
    Ze0      = A * exp( -xi * abs(E0) );
    R1vec2(i)= alpha * Ze0 * ( alpha + (1 - alpha)* k1^(-rho) )^(1/rho - 1);
end;
%% figure
h = figure;
plot(R1vec,mtvec,'-','LineWidth',2,'Color','blue');
hold on;
plot(R1vec2,mtvec2,'--','LineWidth',2,'Color','red');
legend({'Capital Subsidy','Lump-sum'},'Interpreter','latex','FontSize',14,'Location','northeast');
xlabel('Interest Rate $R_{t+1}$','Interpreter','latex','FontSize',14,'Color','black');
ylabel('Mitigation share $m_t$','Interpreter','latex','FontSize',14,'Color','black');
xlim([0 2]);
%dim = [.4 .7 .2 .2];
%str = {'$\theta=$' theta '$\rho =$' rho '$k0 = $' k0 '$A = $' A '$\tau^o_{t+1}=$' tauo1 '$\mathcal{T}_{t+1}=$' T1 };
%annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex');
hold off;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,'R1_mt_plot_cobb.pdf','-dpdf','-r0');