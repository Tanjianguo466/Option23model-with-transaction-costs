function over2_methods()
    % Parameter settings
             
    V0 = 0.25;         % Initial variance
    K = 100;           % Strike price
    r = 0.0;           % Risk-free rate
    kappa = 2.5;       % Variance reversion speed
    theta = 0.16;      % Long-term variance level
    sigma = 0.45;      % Volatility of volatility
    rho = 0.1;         % Correlation between stock and variance
    T = 1;             % Time to maturity
    dt = 1/252;        % Time step (assuming 252 trading days)
    N = 50000;         % Number of samples
    fprintf('Antithetic method (N = %d):\n', N);
    % Antithetic variates method
    for S0 = 80:10:120
       
        [P_antithetic, Var_antithetic] = monte_carlo_3over2_antithetic(S0, V0, K, r, kappa, theta, sigma, rho, T, N, dt);
       
        
        fprintf('  Stock price: %.4f\n', S0);
        fprintf('  Option price: %.4f\n', P_antithetic);
    end
end

function [P, Var] = monte_carlo_3over2_antithetic(S0, V0, K, r, kappa, theta, sigma, rho, T, N, dt)
    % Initialization
    payoff = zeros(N, 1);   % Option payoffs
    steps = round(T / dt);  % Total number of time steps

    % Generate random numbers
    Z1 = randn(N, steps);   % Random terms for stock price
    Z2 = randn(N, steps);   % Random terms for variance

    % Simulate paths
    for i = 1:N
        % Original path
        S = S0;
        V = V0;
        
        % Antithetic path (only reverse Z1)
        S_antithetic = S0;
        V_antithetic = V0;
        
        for t = 1:steps
            % Random terms for original path
            dW1 = sqrt(dt) * Z1(i, t);
            dW2 = sqrt(dt) * (rho * Z1(i, t) + sqrt(1 - rho^2) * Z2(i, t));
            
            % Random terms for antithetic path (only reverse dW1)
            dW1_antithetic = -dW1;
            dW2_antithetic = dW2;  
            
            % Update variance (3/2 model)
            V = V + kappa * (theta - V) * V * dt + sigma * V^(3/2) * dW2;
            V = max(V, 0);
            
            V_antithetic = V_antithetic + kappa * (theta - V_antithetic) * V_antithetic * dt + ...
                          sigma * V_antithetic^(3/2) * dW2_antithetic;
            V_antithetic = max(V_antithetic, 0);
            
            % Update stock price
            S = S + r * S * dt + sqrt(V) * S * dW1;
            S_antithetic = S_antithetic + r * S_antithetic * dt + sqrt(V_antithetic) * S_antithetic * dW1_antithetic;
        end
        
        % Calculate option payoff (average of two paths)
        payoff(i) = 0.5 * (max(S - K, 0) + max(S_antithetic - K, 0));
    end

    % Calculate option price and variance
    P = exp(-r * T) * mean(payoff);
    Var = var(payoff) / N;
end