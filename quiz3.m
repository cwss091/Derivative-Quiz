close all;
clear all;
clc;

mu = 0.05;
sigma = 0.25;
rf_rate = 0.03;
Value_Firm = 100;
Value_Debt = 50;
Num_Stock = 1;
TOM = 2;
T = 1;

%N_Simulation = 10000;
N_Simulation = 10;
dt = 1/52;
Time = 0:dt:T;
T2Maturity = TOM - Time;

K = Value_Firm - Value_Debt; %firm value per - debt value per

dW = sqrt(dt) * randn(length(Time) - 1, N_Simulation);

% S is the firm's value, not the stock price
Log_S_0 = log(Value_Firm) * ones(1, N_Simulation);
d_Log_S_t = (mu-0.5*sigma*sigma)*dt*ones(length(Time)-1, N_Simulation)...
    + sigma * dW;
Log_S_t = [Log_S_0; d_Log_S_t];
Log_S_t = cumsum(Log_S_t);
S_t = exp(Log_S_t);
%1/N_Simulation * sum(S_t(end,:))

[StockPrice,BondPrice] = blsprice(S_t, K, rf_rate, T2Maturity' * ones(1, N_Simulation), sigma, 0);
[StockDelta,BondDelta] = blsdelta(S_t, K, rf_rate, T2Maturity' * ones(1, N_Simulation), sigma, 0);


Stock_P = 1/N_Simulation * sum(StockPrice(1,:))
Bond_P = 50 * exp(-rf_rate * 2) - 1/N_Simulation * sum(BondPrice(1,:))

Stock_D = 1/N_Simulation * sum(StockDelta(1,:))
Bond_D = -1/N_Simulation * sum(BondDelta(1,:))

%% 3


ValueHedged = nan(size(S_t,1), size(S_t,2));
ValueHedged(1, :) = 1000000 * ones(1, length(N_Simulation));

figure;
for i = 1:N_Simulation
    for j = 1:length(Time) - 1
        %dynamic delta hedging
        A = [StockDelta(j, i) -BondDelta(j, i); StockPrice(j, i) 50*exp(-rf_rate*T2Maturity(j))-BondPrice(j, i)];
        b = [0; ValueHedged(j, i)];
        Portf = A \ b;
        ValueHedged(j+1, i) = Portf(1) * StockPrice(j+1, i) + Portf(2) * 50*exp(-rf_rate*T2Maturity(j+1))-BondPrice(j+1, i);
    end
    plot(Time, ValueHedged(:, i));
    hold on;
end

xlabel('Time');
ylabel('Value of Delta Hedged Portfolio');
plot(Time, 1000000*exp(rf_rate*Time),'r')
hold off;


title('Weekly Rebalancing')

    





