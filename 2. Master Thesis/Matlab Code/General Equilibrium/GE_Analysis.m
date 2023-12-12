
clear; clc; close all; 

global gamma beta agrid zgrid gamma11 gamma12 gamma21 gamma22 tol scale_phi
global I J aa zz P


%% Parametrization 

% households 
dr = 0.04;              % discount rate
beta = 1/(1 + dr);      % discount factor 
gamma = 2;              % CRRA coefficient inverse of IES   
phi = -0;                % borrowing limit 
sigma = 0.2;            % variance AR(1) income innovation
rho = 0.9;              % income autoregressive coefficient
scale_phi = 1;           % Scale parameter for the Gumbel Distribution

tol = 1e-13;

% firms 
delta = 0.03;          % depreciation rate 
alpha = 0.34;          % Cobb-douglas capital income share 
omega = 0.33;
eta = 0.33;


% grid size 
I = 1000;                % asset grid 
J = 2;                 % income grid 

gamma11 = 1;
gamma12 = 0.5;
gamma21 = 0.5   ;
gamma22 = 1;

M1 = 0.5;
M2 = 0.5;

%% Grids and discretized income process 

% uniform asset grid 
amax = 5; 
amin = -phi;
agrid = zeros(I,1); 
for i = 1:I 
agrid(i) = amin+(i-1)*(amax-amin)/(I-1);
end 

% for now just use a two point markov process to approximate the log AR(1) process 
% log z' = \rho ln z + \sigma \nu  where \nu \sim N(0,1)

% income grid 
sigmaly2 = (sigma^2)/(1 - rho^2); 
den = exp(1 - sqrt(sigmaly2)) + exp(1 + sqrt(sigmaly2)); 
zl = (2*exp(1-sqrt(sigmaly2)))/den; 
zh = (2*exp(1+sqrt(sigmaly2)))/den;
zgrid = [zl; zh];

% combined grids 
aa = repmat(agrid,1,J); 
zz = repmat(zgrid',I,1);

% transition matrix 
Phi = 1 + log(2/den);
kappa = ((sigma^2)*(1+rho) - (Phi^2)*(1-rho))/(2*(sigma^2)); 
P = [kappa 1-kappa; 1-kappa kappa]; 

% initial guess  
p = [.6; .4]; 
test = 1;
% stationary distribution 
while test > 1e-16    
p1 = P'*p;
test = max(abs(p1 - p));
p = p1;        
end
save('Decision_Rules/agrid.mat','agrid')
%% Steady State General equilibrium

% Initialization 
eps = 1e-02;
it = 0;
itmax = 30; 

% Initial guess on stationary distribution, uniformly distributed
faz11 = ones(I, J)*(1/(I*J));
faz12 = ones(I, J)*(1/(I*J));
faz21 = ones(I, J)*(1/(I*J));
faz22 = ones(I, J)*(1/(I*J));

% Initial guess on steady state occupational probabilities
p11_ss = ones(I, J)*0.5;
p12_ss = ones(I, J)*0.5;
p21_ss = ones(I, J)*0.5;
p22_ss = ones(I, J)*0.5;

% Initial guess on steady state interest rate
rmin = - delta + 1e-8; 
rmax = dr - 1e-8;
r = .5*(rmin+rmax);

tic
while it < itmax

% Here the interest rate is given and we have a guess of the 
% Labor supply. We 
L1 = M1*sum((p11_ss.*faz11.*gamma11), 'all') + M2*sum((p21_ss.*faz21.*gamma21),'all');
L2 = M1*sum((p12_ss.*faz12.*gamma12), 'all') + M2*sum( (p22_ss.*faz22.*gamma22),'all');

% equilibrium prices   
K = (((r + delta)*(L1^(-omega))*(L2^(-eta)))/alpha)^(1/(alpha-1));  
w1 = (omega)*L1^(omega - 1)*(K^alpha)*(L2^eta);
w2 = eta*(K^alpha)*(L1^omega)*(L2^(eta - 1));

% household optimization
[v11, v12, v21, v22, ap11, ap12, ap21, ap22, c11, c12, c21, c22] = VFI_GE(r, w1, w2);

% stationary distribution  
[faz11] = CK_GE(ap11);
[faz12] = CK_GE(ap12);
[faz21] = CK_GE(ap21);
[faz22] = CK_GE(ap22);

p11_ss = exp(v11./scale_phi)./(exp(v11./scale_phi) + exp(v12./scale_phi));
p12_ss = 1 - p11_ss;
p21_ss = exp(v21./scale_phi)./(exp(v21./scale_phi) + exp(v22./scale_phi));
p22_ss = 1 - p21_ss;

% market clearing
K1 = M1*sum(faz11.*p11_ss.*aa + faz12.*p12_ss.*aa,'all') + M2*sum(p21_ss.*faz21.*aa + p22_ss.*faz22.*aa,'all');

res = (K1 - K);

% update r using bisection 
if res > eps
    fprintf('Excess Supply, r = %.4f, w1 = %.4f, w2 = %.4f \n',r,w1, w2);
    rmax = r; 
    r = 0.5*(r+rmin);
    it = it + 1;
elseif res < -eps
    fprintf('Excess Demand, r = %.4f, w1 = %.4f, w2 = %.4f \n',r,w1,w2);
    rmin = r; 
    r = 0.5*(r+rmax);
    it = it + 1;
elseif res < eps
    fprintf('Equilibrium Found, (r, w1, w2, it) = (%.4f,%.4f,%.4f,%.4f) \n',r,w1,w2,it);
    disp('');
    break
end

disp(res); 

end
toc

%% Save Steady State Prices
prices = [r, w1, w2];
save('Decision_Rules\prices.mat','prices')

%% Recompute The equilibrium
tic
[v11, v12, v21, v22, ap11, ap12, ap21, ap22, c11, c12, c21, c22] = VFI_GE(r, w1, w2);
toc

tic
% stationary distribution  
[faz11] = CK_GE(ap11);
[faz12] = CK_GE(ap12);
[faz21] = CK_GE(ap21);
[faz22] = CK_GE(ap22);
toc

% Occupational Probabilities
p11 = exp(v11./scale_phi)./(exp(v11./scale_phi) +exp(v12./scale_phi));
p12 = exp(v12./scale_phi)./(exp(v11./scale_phi) +exp(v12./scale_phi));
p21 = exp(v21./scale_phi)./(exp(v21./scale_phi) + exp(v22./scale_phi));
p22 = exp(v22./scale_phi)./(exp(v21./scale_phi) + exp(v22./scale_phi));

%% Save the equilibrium objects
save('Decision_Rules\v11.mat','v11')
save('Decision_Rules\v12.mat','v12')
save('Decision_Rules\v21.mat','v21')
save('Decision_Rules\v22.mat','v22')

save('Decision_Rules\ap11.mat','ap11')
save('Decision_Rules\ap12.mat','ap12')
save('Decision_Rules\ap21.mat','ap21')
save('Decision_Rules\ap22.mat','ap22')

save('Decision_Rules\c11.mat','c11')
save('Decision_Rules\c12.mat','c12')
save('Decision_Rules\c21.mat','c21')
save('Decision_Rules\c22.mat','c22')

save('Decision_Rules\faz11.mat','faz11')
save('Decision_Rules\faz12.mat','faz12')
save('Decision_Rules\faz21.mat','faz21')
save('Decision_Rules\faz22.mat','faz22')

save('Decision_Rules\p11.mat','p11')
save('Decision_Rules\p12.mat','p12')
save('Decision_Rules\p21.mat','p21')
save('Decision_Rules\p22.mat','p22')

%% Plot Section
figure(1)
title('Value Functions - Type 1')
hold on
plot(agrid, v11(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, v11(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, v12(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, v12(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
hold off
legend

figure(2)
title('Value Functions - Type 2')
hold on
plot(agrid, v21(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, v21(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, v22(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, v22(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
hold off
legend

figure(3)
title('Consumption Functions - Type 1')
hold on
plot(agrid, c11(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, c11(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, c12(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, c12(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
hold off
legend

figure(4)
title('Consumption Functions - Type 2')
hold on
plot(agrid, c21(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, c21(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, c22(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, c22(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
hold off
legend

figure(5)
title('Asset Functions - Type 1')
hold on
plot(agrid, ap11(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, ap11(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, ap12(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, ap12(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
hold off
legend

figure(6)
title('Asset Functions - Type 2')
hold on
plot(agrid, ap21(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, ap21(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, ap22(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, ap22(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
hold off
legend

figure(7)
title('Occupational Probabilities - Type 1')
hold on
plot(agrid, p11(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, p11(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, p12(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, p12(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
hold off
legend

figure(8)
title('Occupational Probabilities - Type 2')
hold on
plot(agrid, p21(:,1),'DisplayName','Occ1 - zl','LineWidth',2,'color',"#ff9696")
plot(agrid, p21(:,2),'DisplayName','Occ1 - zh','LineWidth',2,'color',"#FF0000")
plot(agrid, p22(:,1),'DisplayName','Occ2 - zl','LineWidth',2,'color',"#8293ff")
plot(agrid, p22(:,2),'DisplayName','Occ2 - zh','LineWidth',2,'color',"#0023ff")
legend
hold off


figure(9)
title('Occupation 1')
hold on
fa = sum(faz11,2);
abar = linspace(amin,amax,1000); 
fabar = interp1(agrid,fa,abar); 
b = bar(abar,fabar,'hist'); 
%ylim([0,0.05]);
%xlim([-0.2,agrid(I)]);
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('Wealth distribution','Interpreter','latex','FontSize',14); 
hold off

figure(10)
title('Occupation 2')
hold on
fa = sum(faz12,2);
abar = linspace(amin,amax,1000); 
fabar = interp1(agrid,fa,abar); 
b = bar(abar,fabar,'hist'); 
%ylim([0,0.05]);
%xlim([-0.2,agrid(I)]);
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('Wealth distribution','Interpreter','latex','FontSize',14); 
hold off

figure(11)
title('Total Wealth Distribution')
hold on
fa = sum( (faz11.*p11 + faz12.*p12)*M1 + (faz21.*p21 + faz22.*p22)*M2, 2);
abar = linspace(amin,amax,1000); 
fabar = interp1(agrid,fa,abar); 
b = bar(abar,fabar,'hist'); 
%ylim([0,0.05]);
%xlim([-0.2,agrid(I)]);
xlabel('$a$','Interpreter','latex','FontSize',14); 
ylabel('Wealth distribution','Interpreter','latex','FontSize',14); 
hold off


