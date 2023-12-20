%% Numerical Option Pricing - Antonio Cozzolino P09000193

% Constants
r = 0.02;
mu = 0.15;
sigma = 0.2;
S0 = 100;


%% Exercise 2.1
T = 1;
K = 130;
N = [10,100,1000,10000,100000];

% Since we know the distribution of the stock path is lognormal we can
% compute P(S(T)<130) by using the logncdf function as a theoretical reference.
% The code is divided in two sections: the first one is related to the 
% computation under risk neutral probability, the second one under
% subjective probability. 

% ################################ Risk Neutral Probability ################
% The first variable contain the theoretical probability of the stock price
% to be <K, the second is a vector of 0 that will store the different
% estimated probabilities after each iteration based on the number of
% simualtions.
th_prob_rn = logncdf(K, log(S0) + (r - (sigma^2)/2)*T, sigma*sqrt(T));
est_prob_rn = zeros(1,5);

% The first variable contains the theoretical expected return while the
% second, again, is a vector of 0 that will store the different estimated
% returns after each iteration of the loop.
th_exp_return_rn = S0*exp(r*T)/S0 - 1;
est_exp_return_rn = zeros(1,5);


% The loop iterate on the dimension of N which is the vector of the
% simulations number, numSim sets the number of simulations, then the code
% generates a vector of standard normal distribution for the stochastic
% part of the price. Log price is computed and then at each iteration the
% expected return is computed and put into the vector for the plot, the
% same is done for the probability of the price to be <130, computed as
% number of positive cases/total cases.
for i = 1:length(N)
    numSim = N(i);                     
    Z = normrnd(0,1,numSim,1);         
    logS = log(S0) + (r - (sigma^2)/2)*T  + (sigma*Z*sqrt(T));
    est_exp_return_rn(i) = (mean(exp(logS))-S0)/S0;
    est_prob_rn(i) = sum(exp(logS)<130)/numSim;
end

% This is the PLOT of the Expeted Return
figure(1)
subplot(2,2,1)
plot(1:length(N),est_exp_return_rn,1:length(N),[th_exp_return_rn,th_exp_return_rn,th_exp_return_rn,th_exp_return_rn,th_exp_return_rn])
title('Expected Return - Risk Neutral Probability')
legend('Theoretical','True Value')
xlabel('Number of Simulations')
xticklabels(N)

% This is the PLOT of the Probability of S<K
subplot(2,2,2)
plot(1:length(N),est_prob_rn,1:length(N),[th_prob_rn,th_prob_rn,th_prob_rn,th_prob_rn,th_prob_rn])
title('Probability - Risk Neutral Probability')
legend('Theoretical','True Value')
xlabel('Number of Simulations')
xticklabels(N)

% ################################ Subjective Probability #################
% The first variable contain the theoretical probability of the stock price
% to be <K, the second is a vector of 0 that will store the different
% estimated probabilities after each iteration based on the number of
% simualtions.
theoretical_prob_sbj = logncdf(K, log(S0) + (mu - (sigma^2)/2)*T, sigma*sqrt(T));
estimated_prob_sbj = zeros(1,5);

% The first variable contains the theoretical expected return while the
% second, again, is a vector of 0 that will store the different estimated
% returns after each iteration of the loop.
theoretical_exp_return_sbj = S0*exp(mu*T)/S0 - 1;
estimated_exp_return_sbj = zeros(1,5);

clear logS Z
for i = 1:length(N)
    numSim = N(i);                     
    Z = normrnd(0,1,numSim,1);         
    logS = log(S0) + (mu - (sigma^2)/2)*T  + (sigma*Z*sqrt(T));
    estimated_exp_return_sbj(i) = mean(exp(logS)- S0)/S0;
    estimated_prob_sbj(i) = sum(exp(logS)<130)/numSim;
end

% Expected Return PLOT
subplot(2,2,3)
plot(1:length(N),estimated_exp_return_sbj,1:length(N),[theoretical_exp_return_sbj,theoretical_exp_return_sbj,theoretical_exp_return_sbj,theoretical_exp_return_sbj,theoretical_exp_return_sbj])
title('Expected Return - Subjective Probability')
legend('Theoretical','True Value')
xlabel('Number of Simulations')
xticklabels(N)

% Probability PLOT
subplot(2,2,4)
plot(1:length(N),estimated_prob_sbj,1:length(N),[theoretical_prob_sbj,theoretical_prob_sbj,theoretical_prob_sbj,theoretical_prob_sbj,theoretical_prob_sbj])
title('Probability - Subjective Probability')
legend('Theoretical','True Value')
xlabel('Number of Simulations')
xticklabels(N)

%% Exercise 2.2: Compute numerically the price of a "down-and-out" put option with maturity T.

% The payoff of the option is the following, it is worthless if at any time
% before maturity the asset price falls below the level S_b, otherwise the
% final payoff is max(K-S_T,0). So basically the final payoff will be the
% usual max(K-S_T,0) conditional on the fact that overtime, before maturity
% the stock price never goes below 0.95

S0 = 1;
K = 1;
S_b = 0.95;
dt = 1/250;
days = 250;


% Here we set the number of simulation to be 1000 for speed reason, higher
% number of simulations leads to higher precision, however this is
% computationally hard since every stock path must be created and then we
% access to the array checking all the conditions.
numSim = 1000;
stock_path = zeros(numSim,days);
stock_path(:,1) = S0;
    for i=1:numSim-1
        for j=1:days-1
            stock_path(i,j+1)=stock_path(i,j)*(1+ r*dt)+sigma*stock_path(i,j)*sqrt(dt)*normrnd(0,1,1,1);
        end
    end 

% Plot of stock path for check that the simulation of the trajectory is ok
% figure(2)
% for i=1:numSim-1
%     plot(1:days,stock_path(i,:),1:days,repelem(S_b,days),'red')
%     hold on
% end
% title('Trajectory plot')
% ylabel('Price')
% xlabel('Time')

% At this point we instantiate the payoff vector, the algorithm works as
% follows: it checks into the stock_path matrix if any entry if <0.95, if
% so it jumps all the computation and assign 0 as payoff (for speed reason,
% once we reached the point for which the option expire it is useless to go
% further into the stock path), otherwise it goes on until the terminal 
% value where max(K-S(T),0) is computed. 
% Then, after we have the payoff vector, we compute the mean and discount 
% for get the final put price.

payoff_vector_put = zeros(1,numSim);
for i=1:numSim-1
    for j=1:days-1
        if stock_path(i,j) < S_b
            payoff_vector_put(i) = 0;
            continue
        elseif stock_path(i,j) >= S_b
            if j==days-1
                payoff_vector_put(i) = max( K - stock_path(i,j),0);
            end
        end
    end
end

% Here is the final price of the put printed at display.
pricePUT = exp(-r*T)*mean(payoff_vector_put);
disp('Price of the put:');
disp(pricePUT)


