clear;
clc;

% Parameters
alpha    = 1/3;
beta     = .7;
delta    = 0;
gamma    = 1;
xi       = 1;
m0       = 0;
m1       = 0;
tauo0    = 0;
tauo1    = 0;
tauo2    = 0;
k0       = 2;
M0       = 0;
Emin     = 0;
zeta     = 1;
A        = 3;

% Number of iterations and stepsize
N  = 56;
di = 0.01;

% Initial guess
guess0 = 0;
guess1 = 0;

I0  = A * exp(-xi * Emin) * ( 1 - alpha ) * k0^alpha * ( 1 - alpha / ( 1 - alpha ) * tauo0 );
k1  = beta / ( 1 + beta ) * ( 1 - m0) * I0;
E0  = ( 1 - delta ) * Emin + zeta*k0 - gamma * M0;
R1  = alpha * A * exp(-xi * E0) * k1^(alpha-1);

for i =  1:N  

parameters = [ alpha beta delta gamma xi m0+(i-1)*di m1 tauo0 tauo1 tauo2 k0 M0 Emin zeta A];

[tauo10vec(i), fval, eflag, out] = fsolve(@(XX) dW0(XX, parameters),[guess0],...
optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000));

[tauo11vec(i), fval, eflag, out] = fsolve(@(XX) dW1(XX, parameters),[guess1],...
optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,'MaxIter', 100000));

guess0    = tauo10vec(i);
guess1    = tauo11vec(i);
mt(i)   = m0+(i-1)*di;

parameters1 = [ alpha beta delta gamma xi m0+(i-1)*di m1 tauo0 tauo10vec(i) tauo2 k0 M0 Emin zeta A];
parameters2 = [ alpha beta delta gamma xi m0+(i-1)*di m1 tauo0 tauo11vec(i) tauo2 k0 M0 Emin zeta A];

[U00(i), U10(i), k10(i), k20(i)] = welfarelog(parameters1);
[U01(i), U11(i), k11(i), k21(i)] = welfarelog(parameters2);

end

% Figures

h = figure;
plot(mt,tauo10vec,'-','LineWidth',2,'Color','blue');
hold on;
plot(mt,tauo11vec,'--','LineWidth',2,'Color','red');
[xout,yout] = intersections(mt(2:N),tauo10vec(2:N),mt(2:N),tauo11vec(2:N),1);
mtmatrix    = [mt ; tauo10vec ; tauo11vec];
mtout       = mtmatrix(:,mtmatrix(1,:)<=xout);
Tout10      = mtout(2,:);
Tout11      = mtout(3,:);
x2          = [mtout(1,:), fliplr(mtout(1,:))];
inBetween   = [Tout10, fliplr(Tout11)];
pareto      = fill(x2, inBetween, 'g','LineStyle','none');
set(pareto,'facealpha',.15);
legend({'$\Omega^c$','$\psi^c$'},'Interpreter','latex','FontSize',18,'Location','northwest');
grid on;
hold off;
% hold on;
% hold off;
xlabel('Mitigation share \hspace{0.15cm} $m_t$','Interpreter','latex','FontSize',14,'Color','black');
ylabel('Transfer rate \hspace{0.15cm} $\tau^o_{t+1}$','Interpreter','latex','FontSize',14,'Color','black');
%legend({'FT','ETS2'},'Location','northwest','Orientation','horizontal')
xlim([0 0.55]);
ylim([0 1]);
xticks(0:0.05:0.55);
yticks(0:0.25:1.5);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,'Cobb_Subsidy_R1_lumsuper','-dpdf','-r0');