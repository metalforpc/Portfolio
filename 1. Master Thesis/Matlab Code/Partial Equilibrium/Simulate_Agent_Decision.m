%% Load variables
load('Decision_Rules/ap11.mat')
load('Decision_Rules/ap12.mat')
load('Decision_Rules/c11.mat')
load('Decision_Rules/c12.mat')
load('Decision_Rules/faz11.mat')
load('Decision_Rules/faz12.mat')
load('Decision_Rules/p11.mat')
load('Decision_Rules/p12.mat')
load('Decision_Rules/v11.mat')
load('Decision_Rules/v12.mat')
load('Decision_Rules/agrid.mat')

%% Simulate Agent Decision Path
% This code simulates a trajectory for an individual subject to 
% taste shock and idiosyncratic shock.

% We want to record the following variables

% taste shock, idiosyncratic shock, asset, consumption, value,
% occupational_choice 

T = 100; % Number of time periods
vals = zeros(T, 6); % values of the simulation, Tx6 array

