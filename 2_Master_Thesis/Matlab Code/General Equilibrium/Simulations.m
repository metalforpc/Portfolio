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

% Deccomment for taste shocks
%phi_mat = evrnd(0, scale_phi, T+1, (N1+N2)*2);

% No taste shock
phi_mat = zeros(T+1, (N1+N2)*2);
z_mat = zeros(T+1, N1+N2);

% Simulate the time series for every individual 
for i=1:N1+N2
    z_mat(:,i) = simulate(mc, T);
end

% ALTERNATIVE: DE-COMMENT to obtain a series with no shocks
% phi_mat = zeros(T+1,(N1+N2)*2);
% z_mat = ones(T+1,N1+N2);

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
C = sum(consumption_mat,2)/(N1+N2);
K = sum(asset_mat,2)/(N1+N2);
m11 = sum(occ_mat(:,1:N1) == 1,2)/N1;
m12 = sum(occ_mat(:,1:N1) == 2,2)/N1;
m21 = sum(occ_mat(:,N1+1:N1+N2) == 1,2)/N2;
m22 = sum(occ_mat(:,N1+1:N1+N2) == 2,2)/N2;

m1 = M1*m11 + M2*m21;
m2 = M1*m12 + M2*m22;

L1 = gamma11*m11*M1 + gamma21*m21*M2;
L2 = gamma12*m12*M1 + gamma22*m22*M2;

e1 = m11*M1 + m21*M2;
e2 = m12*M1 + m22*M2;

gammabar1 = L1./e1;
gammabar2 = L2./e2;

L = L1 + L2;

Y = (K.^alpha).*(L1.^(omega)).*(L2.^(eta));

I = Y - C;

%% Prices 
r = alpha.*(K.^(1-alpha)).*(L1.^(omega)).*(L2.^(eta));
w1 = omega.*(K.^(alpha)).*(L1.^(1 - omega)).*(L2.^(eta));
w2 = eta.*(K.^(alpha)).*(L1.^(omega)).*(L2.^(1 -  eta));
w = w1.*gammabar1.*e1 + w2.*gammabar2.*e2;

% %% HP DECOMPOSITION AND CORRELATIONS
% [trendC, cycleC] = hpfilter(C(1:T-1));
% [trendK, cycleK] = hpfilter(K(1:T-1));
% [trendm11, cyclem11] = hpfilter(m11(1:T-1));
% [trendm12, cyclem12] = hpfilter(m12(1:T-1));
% [trendm21, cyclem21] = hpfilter(m21(1:T-1));
% [trendm22, cyclem22] = hpfilter(m22(1:T-1));
% [trendL1, cycleL1] = hpfilter(L1(1:T-1));
% [trendL2, cycleL2] = hpfilter(L2(1:T-1));
% [trendr, cycler] = hpfilter(r(1:T-1));
% [trendw1, cyclew1] = hpfilter(w1(1:T-1));
% [trendw2, cyclew2] = hpfilter(w2(1:T-1));
% [trendY, cycleY] = hpfilter(Y(1:T-1));
% [trendI, cycleI] = hpfilter(I(1:T-1));
% [trendL, cycleL] = hpfilter(L(1:T-1));
% [trendw, cyclew] = hpfilter(w(1:T-1));
% [trende1, cyclee1] = hpfilter(e1(1:T-1));
% [trende2, cyclee2] = hpfilter(e1(1:T-1));
% 
% %% Data matrix to compute correlation on
% cormat = [cycleC, cycleK, cyclem11, cyclem12, cyclem21, cyclem22, cycleL1, cycleL2, cycler, cyclew1, cyclew2, cycleY, cycleI, cycleL, cyclew];
% 
% correlations = corrcoef(cormat);
% lbls = {'C', 'K', 'm_{11}', 'm_{12}', 'm_{21}', 'm_{22}', 'L_1', 'L_2', 'r', 'w_1', 'w_2', 'Y', 'I', 'L', 'w'};
% lbls = categorical(lbls);
% heatmap(correlations, 'XDisplayLabels', lbls, 'YDisplayLabels', lbls, 'Colormap', autumn)
% 
% %% Distributional Properties 
% 
% figure(1)
% hold on 
% ax = gca;
% set(gca,'TickLabelInterpreter', 'latex');
%     set(gca,...
%             'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%          'FontName','cmr10',...
%         'FontSize',14,...
%         'Box','off'); 
% 
% subplot(2,3,1)
% plot(K(1:200))
% xlabel('$t$','Interpreter','latex')
% ylabel('$K_t$','Interpreter','latex')
% 
% subplot(2,3,2)
% plot(m11(1:200))
% xlabel('$t$','Interpreter','latex')
% ylabel('${M_{1,1}}_t$','Interpreter','latex')
% 
% subplot(2,3,3)
% plot(m12(1:200))
% xlabel('$t$','Interpreter','latex')
% ylabel('${M_{1,2}}_t$','Interpreter','latex')
% 
% subplot(2,3,4)
% histogram(asset_mat(2,:), 50, 'Normalization', 'pdf')
% ylabel('$t=2$','Interpreter','latex')
% 
% subplot(2,3,5)
% histogram(asset_mat(10,:), 50, 'Normalization', 'pdf')
% ylabel('$t=10$','Interpreter','latex')
% 
% subplot(2,3,6)
% histogram(asset_mat(T,:), 50, 'Normalization', 'pdf')
% ylabel('$t=10000$','Interpreter','latex')
% hold off
% 
% 
% % 
% % for i=1:T
% %     histogram(asset_mat(i,:), 10, 'Normalization', 'pdf')
% %     pause(0.1)
% % end
% 
% %% Correlations
% corr_shock_occ = [];
% pval = [];
% 
% for i=1:N1+N2
%     [corr_shock_occ(i), pval(i)] = corr(asset_mat(:,i), occ_mat(:,i));
% end
% 
% %%
% 
% [co, pval] = corr((cyclee1+cyclee2), cyclew)