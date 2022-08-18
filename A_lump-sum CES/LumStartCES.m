%%
clear;
clc;

% Parameters
alpha    = 1/3;     %parameters(1)
beta     = .7;     %parameters(2)
delta    = 0;       %parameters(3)
gamma    = 1;       %parameters(4)
rho      = 0.0001;     %parameters(5)
theta    = 0.01;     %parameters(6)
A        = 3;       %parameters(7)
m0       = 0;       %parameters(8) No social contract case
m1       = 0;       %parameters(9)
T0       = 0;       %parameters(10) The lump sum transfer
T1       = 0;       %parameters(11) No social contract case
T2       = 0;       %parameters(12)
k0       = 2;       %parameters(13)
M0       = 0;       %parameters(14)
Emin     = 0;       %parameters(15)
zeta     = 1;       %parameters(16)
xi       = 1;       %parameters(17)

parameters = [ alpha beta delta gamma rho theta A m0 m1 T0 T1 T2 k0...
    M0 Emin zeta xi];

N        = 56; % Number iteration for social contract
di       = 0.01; % stepsize for social contract m

k1 = fsolve(@(X) lumsolvek1(X, parameters), 5);
Zemin    = A * exp(-xi * abs(Emin) );
w0       = (1 - alpha) * Zemin * ( alpha * k0^rho + (1 - alpha) )^(1/rho - 1);
E0       = (1 - delta) * Emin + zeta*k0 - gamma * M0;
Ze0      = A * exp( -xi * abs(E0) );
R1       = alpha * Ze0 * ( alpha + (1 - alpha)* k1^(-rho) )^(1/rho - 1);

%% 
 T10vec = zeros(1,N); %Preallocation for Lump Sum T1 generation 0
 T11vec = zeros(1,N); %Preallocation for Lump Sum T1 generation 1
 mt     = zeros(1,N); %Preallocation for mitigation
 guessT10 = 0.1;
 guessT11 = 0.1;

for i = 1:N
    m0s       =(i-1)*di; %the simulated contract for mitigation
    [T10vec(i), fval, eflag, out] = fsolve(@(x) lumsolvedW0(x, parameters, m0s),...
        guessT10,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    guessT10 = T10vec(i);
    [T11vec(i), fval, eflag, out] = fsolve(@(x) lumsolvedW1(x, parameters, m0s),...
        guessT11,optimset('TolX',1.0e-8,'TolFun',1.0e-8,'MaxFunEvals',100000,...
        'MaxIter', 100000));
    guessT11 = T11vec(i);
    mt(i)  = (i-1)*di;
end;
%% figure
h = figure;
plot(mt,T10vec,'-','LineWidth',2,'Color','blue');
hold on;
plot(mt, T11vec,'--','LineWidth',2,'Color','red');
[xout,yout] = intersections(mt(2:N),T10vec(2:N),mt(2:N),T11vec(2:N),1);
mtmatrix    = [mt ; T10vec ; T11vec];
mtout       = mtmatrix(:,mtmatrix(1,:)<=xout);
Tout10      = mtout(2,:);
Tout11      = mtout(3,:);
x2          = [mtout(1,:), fliplr(mtout(1,:))];
inBetween   = [Tout10, fliplr(Tout11)];
pareto      = fill(x2, inBetween, 'g','LineStyle','none');
set(pareto,'facealpha',.15);
grid on;
%xticks(0:0.0025:0.025);
%yticks(0:0.02:0.2);
ylim([0 inf]);
xlim([0 0.55]);
legend({'$\Omega^l$','$\psi^l$'},'Interpreter','latex','FontSize',18,'Location','northwest');
xlabel('Mitigation share \hspace{0.15cm}$m_t$','Interpreter','latex','FontSize',14,'Color','black');
ylabel('Lump-sum transfer \hspace{0.15cm} $\mathcal{T}_{t+1}$','Interpreter','latex','FontSize',14,'Color','black');
hold off;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,'11LumCES_k01_A8_thetamin1_highR1_lowPOL001.pdf','-dpdf','-r0');