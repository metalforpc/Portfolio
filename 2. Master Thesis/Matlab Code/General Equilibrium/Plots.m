%% Load Functions
load('Decision_Rules\v11.mat','v11')
load('Decision_Rules\v12.mat','v12')
load('Decision_Rules\v21.mat','v21')
load('Decision_Rules\v22.mat','v22')

load('Decision_Rules\ap11.mat','ap11')
load('Decision_Rules\ap12.mat','ap12')
load('Decision_Rules\ap21.mat','ap21')
load('Decision_Rules\ap22.mat','ap22')

load('Decision_Rules\c11.mat','c11')
load('Decision_Rules\c12.mat','c12')
load('Decision_Rules\c21.mat','c21')
load('Decision_Rules\c22.mat','c22')

load('Decision_Rules\faz11.mat','faz11')
load('Decision_Rules\faz12.mat','faz12')
load('Decision_Rules\faz21.mat','faz21')
load('Decision_Rules\faz22.mat','faz22')

load('Decision_Rules\p11.mat','p11')
load('Decision_Rules\p12.mat','p12')
load('Decision_Rules\p21.mat','p21')
load('Decision_Rules\p22.mat','p22')

%% Parameters
% households 
dr = 0.04;              % discount rate
beta = 1/(1 + dr);      % discount factor 
gamma = 2;              % CRRA coefficient inverse of IES   
phi = -0;                % borrowing limit 
sigma = 0.2;            % variance AR(1) income innovation
rho = 0.9;              % income autoregressive coefficient
scale_phi = 1;           % Scale parameter for the Gumbel Distribution

tol = 1e-08;

% firms 
delta = 0.03;          % depreciation rate 
alpha = 0.34;          % Cobb-douglas capital income share 
omega = 0.33;
eta = 0.33;

% grid size 
I = 1000;                % asset grid 
J = 2;                 % income grid 

gamma11 = 1.5;
gamma12 = 0.5;
gamma21 = 0.5;
gamma22 = 1.5;

M1 = 0.5;
M2 = 0.5;

% uniform asset grid 
amax = 10; 
amin = -phi;
agrid = zeros(I,1); 
for i = 1:I 
agrid(i) = amin+(i-1)*(amax-amin)/(I-1);
end 

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