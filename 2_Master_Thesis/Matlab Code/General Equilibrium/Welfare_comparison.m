%%
close all 
clear all

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
for s=1:2

    T = 10000; % Number of time periods
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
    
    asset_mat(1,:) = agrid(2); % Same initial condition for everyone
    
    % First we have to draw a time series for the productivity shock and for
    % phi. Phi is iid distributed and have no serial correlation is sufficient
    % to draw a Tx(N1+N2*2) sample from a gumbel distribution. Recall that in
    % this matrix the row is the time index while each couple of column is a
    % preference shock for an individual
    
    % In the first loop we have taste shocks
    if s == 1
        phi_mat = evrnd(0, scale_phi, T+1, (N1+N2)*2);
    end 
    
    % No taste shock
    %phi_mat = zeros(T+1, (N1+N2)*2);
    z_mat = zeros(T+1, N1+N2);
    
    % Simulate the time series for every individual 
    for i=1:N1+N2
        z_mat(:,i) = simulate(mc, T);
    end
    
    % In the second loop no taste shocks
    if s == 2
        phi_mat = zeros(T+1,(N1+N2)*2);
    end

    %z_mat = ones(T+1,N1+N2);
    
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
    if s == 1
        C_1 = sum(consumption_mat,2)/(N1+N2);
        K_1 = sum(asset_mat,2)/(N1+N2);
        m11_1 = sum(occ_mat(:,1:N1) == 1,2)/N1;
        m12_1 = sum(occ_mat(:,1:N1) == 2,2)/N1;
        m21_1 = sum(occ_mat(:,N1+1:N1+N2) == 1,2)/N2;
        m22_1 = sum(occ_mat(:,N1+1:N1+N2) == 2,2)/N2;
        
        m1_1 = M1*m11_1 + M2*m21_1;
        m2_1 = M1*m12_1 + M2*m22_1;
        
        L1_1 = gamma11*m11_1*M1 + gamma21*m21_1*M2;
        L2_1 = gamma12*m12_1*M1 + gamma22*m22_1*M2;
        
        Y_1 = (K_1.^alpha).*(L1_1.^(omega)).*(L2_1.^(eta));
        
        I_1 = Y_1 - C_1;
        
        %% Prices 
        r_1 = alpha.*(K_1.^(1-alpha)).*(L1_1.^(omega)).*(L2_1.^(eta));
        w1_1 = omega.*(K_1.^(alpha)).*(L1_1.^(1 - omega)).*(L2_1.^(eta));
        w2_1 = eta.*(K_1.^(alpha)).*(L1_1.^(omega)).*(L2_1.^(1 -  eta));
    else
        C_2 = sum(consumption_mat,2)/(N1+N2);
        K_2 = sum(asset_mat,2)/(N1+N2);
        m11_2 = sum(occ_mat(:,1:N1) == 1,2)/N1;
        m12_2 = sum(occ_mat(:,1:N1) == 2,2)/N1;
        m21_2 = sum(occ_mat(:,N1+1:N1+N2) == 1,2)/N2;
        m22_2 = sum(occ_mat(:,N1+1:N1+N2) == 2,2)/N2;
        
        m1_2 = M1*m11_2 + M2*m21_2;
        m2_2 = M1*m12_2 + M2*m22_2;
        
        L1_2 = gamma11*m11_2*M1 + gamma21*m21_2*M2;
        L2_2 = gamma12*m12_2*M1 + gamma22*m22_2*M2;
        
        Y_2 = (K_2.^alpha).*(L1_2.^(omega)).*(L2_2.^(eta));
        
        I_2 = Y_2 - C_2;
        
        %% Prices 
        r_2 = alpha.*(K_2.^(1-alpha)).*(L1_2.^(omega)).*(L2_2.^(eta));
        w1_2 = omega.*(K_2.^(alpha)).*(L1_2.^(1 - omega)).*(L2_2.^(eta));
        w2_2 = eta.*(K_2.^(alpha)).*(L1_2.^(omega)).*(L2_2.^(1 -  eta));
    end 
end
%% Plots
figure(1)
hold on
plot(Y_1(1:100),'DisplayName','$Y_{T}$')
plot(Y_2(1:100),'DisplayName','$Y_{NT}$')
xlabel('$t$','Interpreter','latex')
ylabel('$Y_t$','Interpreter','latex')
ax = gca;
set(gca,'TickLabelInterpreter', 'latex');
    set(gca,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

figure(2)
hold on
plot(C_1(1:100),'DisplayName','$C_{T}$')
plot(C_2(1:100),'DisplayName','$C_{NT}$')
xlabel('$t$','Interpreter','latex')
ylabel('$C_t$','Interpreter','latex')
ax = gca;
set(ax,'TickLabelInterpreter', 'latex');
    set(ax,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

figure(3)
hold on
plot(L1_1(1:100),'DisplayName','$L1_{T}$')
plot(L1_2(1:100),'DisplayName','$L2_{NT}$')
ax = gca;
set(ax,'TickLabelInterpreter', 'latex');
    set(ax,...
            'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
         'FontName','cmr10',...
        'FontSize',14,...
        'Box','off');
hold off
hl = legend('show','Location','southeast');
set(hl, 'Interpreter','latex')

%%
