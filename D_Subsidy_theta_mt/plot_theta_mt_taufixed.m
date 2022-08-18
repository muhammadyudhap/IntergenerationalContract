%%
clear;
clc;

% Parameters
alpha    = 1/3;     %parameters(1)
beta     = .7;     %parameters(2)
delta    = 0;       %parameters(3)
gamma    = 1;       %parameters(4)
rho      = 0;     %parameters(5)
theta    = 0;     %parameters(6)
A        = 3;       %parameters(7)
m0       = 0;       %parameters(8) No social contract case
m1       = 0;       %parameters(9)
tauo0    = 0;       %parameters(10)
tauo1    = 0.01;       %parameters(11) No social contract case
tauo2    = 0;       %parameters(12)
k0       = 2;       %parameters(13)
M0       = 0;       %parameters(14)
Emin     = 0;       %parameters(15)
zeta     = 1;     %parameters(16)
xi       = 1;       %parameters(17)

parameters = [ alpha beta delta gamma rho theta A m0 m1 tauo0 tauo1 tauo2 k0...
    M0 Emin zeta xi];

N        = 500; % Number iteration 
N2       = 100;
di       = 0.01; % stepsize 

%% 
 guessmt = 0.1; 
for i = 1:N
    thetasim       = i*(-di); %the simulated contract for mitigation
    [mtvec(i), fval, eflag, out] = fsolve(@(x) solvingdW0(x, parameters, thetasim),...
        guessmt,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    guessmt = mtvec(i);
    thetavec(i)  = i*(-di);
end;
for i = 1:N2
    thetasim2       = i*(di); 
    [mtvec2(i), fval, eflag, out] = fsolve(@(x) solvingdW0(x, parameters, thetasim2),...
        guessmt,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    guessmt     = mtvec2(i);
    thetavec2(i)  = i*(di);
end;
%%
mtvecall = [mtvec, NaN, mtvec2];
%mtvecall(mtvecall==0)= NaN;
thetavecall = [thetavec,0, thetavec2];
%thetavecall(thetavecall==0)= NaN;

%% figure
h = figure;
plot(thetavecall,mtvecall,'-','LineWidth',2,'Color','blue');
hold on;
grid on;
%xticks(0:0.05:0.55);
%yticks(0:0.02:0.2);
%ylim([0 1]);
%xlim([0 0.55]);
xlabel('Intertemporal substitutability $\theta$','Interpreter','latex','FontSize',14,'Color','black');
ylabel('Mitigation share $m_t$','Interpreter','latex','FontSize',14,'Color','black');
hold off;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,'sub_plot_theta_mt_taufix_fail.pdf','-dpdf','-r0');
