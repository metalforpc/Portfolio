% %% Load Functions
% load('Decision_Rules\v11.mat','v11')
% load('Decision_Rules\v12.mat','v12')
% load('Decision_Rules\v21.mat','v21')
% load('Decision_Rules\v22.mat','v22')
% 
% load('Decision_Rules\ap11.mat','ap11')
% load('Decision_Rules\ap12.mat','ap12')
% load('Decision_Rules\ap21.mat','ap21')
% load('Decision_Rules\ap22.mat','ap22')
% 
% load('Decision_Rules\c11.mat','c11')
% load('Decision_Rules\c12.mat','c12')
% load('Decision_Rules\c21.mat','c21')
% load('Decision_Rules\c22.mat','c22')
% 
% load('Decision_Rules\faz11.mat','faz11')
% load('Decision_Rules\faz12.mat','faz12')
% load('Decision_Rules\faz21.mat','faz21')
% load('Decision_Rules\faz22.mat','faz22')
% 
% load('Decision_Rules\p11.mat','p11')
% load('Decision_Rules\p12.mat','p12')
% load('Decision_Rules\p21.mat','p21')
% load('Decision_Rules\p22.mat','p22')

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
amax = 5; 
amin = -phi;
agrid = zeros(I,1); 
for i = 1:I 
agrid(i) = amin+(i-1)*(amax-amin)/(I-1);
end 


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

% Create the object for the markov chain
mc = dtmc(P);

%% Simulation

T = 200; % Number of time periods
N1 = 100; % Number of agent of type 1
N2 = 100; % Number of agent of type 2


% Now, recall that we want to generate a time series for each individuals,
% and we want to track Asset, Consumption, Value, Z, phi, occupation
% choice for each one, therefore we construct separate matrix for each 
% variable where each raw is the time period and each column is the individual
asset_mat = zeros(T+1, N1+N2);
consumption_mat = zeros(T+1, N1+N2);
value_mat = zeros(T+1, N1+N2);
occ_mat = zeros(T+1, N1+N2);

% Start from the steady state per-capita asset 
asset_mat(1,:) = agrid(103); % Same initial condition for everyone

% Low productivity state and 0 preference shock
phi_mat = zeros(T+1,(N1+N2)*2);
z_mat = ones(T+1,N1+N2);

%% PREFERENCE SHOCK
% I want to simulate a preference shock for the 30% of type 1 in occupation
% 2
pct = 0.6;
phi_mat(60, N1+2:2:N1+N2*pct+1) = 2;

%% Computation
% Now that we have our shocks, we have to simulate our time series 
% For each agent type and for each time period we have the following loop

% First compute type 1 
for i=1:N1
    for t=1:T
        % Check the productivity state
        prod_state = z_mat(t,i);
        phi1 = phi_mat(t,i);
        phi2 = phi_mat(t,i+1);
        
        % Get current asset id
        current_ass_idx = find(asset_mat(t,i) == agrid);

        % Get a matrix of value
        v = [v11(current_ass_idx,prod_state) + phi1, v12(current_ass_idx, prod_state) + phi2];
        [V, occ] = max(v);
        value_mat(t,i) = V;
        occ_mat(t,i) = occ;

        if occ == 1
            consumption_mat(t,i) = c11(current_ass_idx,prod_state);
            asset_mat(t+1, i) = ap11(current_ass_idx,prod_state);
        else 
            consumption_mat(t,i) = c12(current_ass_idx,prod_state);
            asset_mat(t+1, i) = ap12(current_ass_idx,prod_state);
        end

        

    end
end

% Then for occ 2
for i=N1+1:N1+N2
    for t=1:T
        % Check the productivity state
        prod_state = z_mat(t,i);
        phi1 = phi_mat(t,i);
        phi2 = phi_mat(t,i+1);
        
        % Get current asset id
        current_ass_idx = find(asset_mat(t,i) == agrid);

        % Get a matrix of value
        v = [v21(current_ass_idx,prod_state) + phi1, v22(current_ass_idx, prod_state) + phi2];
        [V, occ] = max(v);
        value_mat(t,i) = V;
        occ_mat(t,i) = occ;

        if occ == 1
            consumption_mat(t,i) = c21(current_ass_idx,prod_state);
            asset_mat(t+1, i) = ap21(current_ass_idx,prod_state);
        else 
            consumption_mat(t,i) = c22(current_ass_idx,prod_state);
            asset_mat(t+1, i) = ap22(current_ass_idx,prod_state);
        end
    end
end

%% Aggregate Measures
C = sum(consumption_mat, 2)/(N1+N2);
K = sum(asset_mat, 2)/(N1+N2);
m11 = sum(occ_mat(:,1:N1) == 1,2)/N1;
m12 = sum(occ_mat(:,1:N1) == 2,2)/N1;
m21 = sum(occ_mat(:,N1+1:N1+N2) == 1,2)/N2;
m22 = sum(occ_mat(:,N1+1:N1+N2) == 2,2)/N2;

L1 = gamma11*zl*m11 + gamma21*zl*m21;
L2 = gamma12*zl*m12 + gamma22*zl*m22;

Y = (K.^alpha).*(L1.^(omega)).*(L2.^(eta));

%% Prices 
r = alpha.*(K.^(1-alpha)).*(L1.^(omega)).*(L2.^(eta));
w1 = omega.*(K.^(alpha)).*(L1.^(1 - omega)).*(L2.^(eta));
w2 = eta.*(K.^(alpha)).*(L1.^(omega)).*(L2.^(1 -  eta));

%% Deviation from SS measures
m12(m12 == 0) = 0.000001;
m21(m21 == 0) = 0.000001;


Cdev = (C - C(50))/C(50);
Kdev = (K - K(50))/K(50);
rdev = (r - r(50))/r(50);
w1dev = (w1 - w1(50))/w1(50);
w2dev = (w2 - w2(50))/w2(50);

m11dev = (m11 - m11(50))/m11(50);
m12dev = (m12 - m12(50))/m12(50);
m21dev = (m21 - m21(50))/m21(50);
m22dev = (m22 - m22(50))/m22(50);

L1dev = (L1 - L1(50))/L1(50);
L2dev = (L2 - L2(50))/L2(50);

Ydev = (Y - Y(50))/Y(50);

%% Plots
figure(1)
hold on 
title('$C$','Interpreter','latex')
plot(Cdev(56:100))
hold off

figure(2)
title('$K$','Interpreter','latex')
hold on 
plot(Kdev(56:100))
hold off

figure(3)
title('$r$','Interpreter','latex')
hold on 
plot(rdev(56:100))
hold off

figure(4)
title('$w_1$','Interpreter','latex')
hold on 
plot(w1dev(56:100))
hold off

figure(5)
title('$w_2$','Interpreter','latex')
hold on 
plot(w2dev(56:100))
hold off

figure(6)
title('$m_{1,1}$','Interpreter','latex')
hold on 
plot(m11dev(56:100))
hold off

figure(7)
title('$m_{1,2}$','Interpreter','latex')
hold on 
plot(m12dev(56:100))
hold off

figure(8)
title('$m_{2,1}$','Interpreter','latex')
hold on 
plot(m21dev(56:100))
hold off

figure(9)
title('$m_{2,2}$','Interpreter','latex')
hold on 
plot(m22dev(56:100))
hold off

figure(10)
title('$L_1$','Interpreter','latex')
hold on 
plot(L1dev(56:100))
hold off

figure(11)
title('$L_2$','Interpreter','latex')
hold on 
plot(L2dev(56:100))
hold off

figure(12)
title('$Y_t$','Interpreter','latex')
hold on 
plot(Ydev(56:100))
hold off
