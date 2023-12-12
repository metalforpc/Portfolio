clear; clc; close all; 

global gamma beta agrid zgrid gamma11 gamma12 scale_phi
global I J aa zz P


%% Parametrization 

% households 
dr = 0.04;              % discount rate
beta = 1/(1 + dr);      % discount factor 
gamma = 2;              % CRRA coefficient inverse of IES   
phi = -0;                % borrowing limit 
sigma = 0.2;            % variance AR(1) income innovation
rho = 0.9;              % income autoregressive coefficient
scale_phi = 1;          % Scale parameter for the Gumbel Distribution
gamma11 = 1.5;            % Productivity in Occupation 1
gamma12 = 0.5;          % Productivity in Occupation 2

% Prices 
w1 = 1; 
w2 = 1;
r = .04;   

% firms 
delta = 0.03;          % depreciation rate 
alpha = 0.34;          % Cobb-douglas capital income share 
omega = 0.33;
eta = 0.33;

% grid size 
I = 1000;                % asset grid 
J = 2;                 % income grid 


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

%% Partial equilibrium 
tic
% household optimization
[v11, v12, ap11, ap12, c11, c12] = VFI_PE(r,w1, w2);
toc

tic
% stationary distribution  
[faz11] = CK_PE(ap11);
[faz12] = CK_PE(ap12);
toc

% Occupational Probabilities
p11 = exp(v11./scale_phi)./(exp(v11./scale_phi) +exp(v12./scale_phi));
p12 = exp(v12./scale_phi)./(exp(v11./scale_phi) +exp(v12./scale_phi));

%% Save all the functions
% Stored in the folder Decision Rules, I do this because I want to take
% separate the plot file and the simulation file from the computational one

% Value functions
save('Decision_Rules/v11.mat','v11')
save('Decision_Rules/v12.mat','v12')

% Asset functions
save('Decision_Rules/ap11.mat','ap11')
save('Decision_Rules/ap12.mat','ap12')

% Consumption functions
save('Decision_Rules/c11.mat','c11')
save('Decision_Rules/c12.mat','c12')

% Stationary Distributions
save('Decision_Rules/faz11.mat','faz11')
save('Decision_Rules/faz12.mat','faz12')

% Occupational Probabilities 
save('Decision_Rules/p11.mat','p11')
save('Decision_Rules/p12.mat','p12')

% Asset Grid
save('Decision_Rules/agrid.mat','agrid')