function prior = prior_Minnesota_struc(n, p, lambda, sig2)
% Constructs Minnesota prior on structural VAR parameters.
%
% INPUTS:
%   n      - number of variables
%   p      - number of lags
%   lambda - vector of hyperparameters [lambda1, lambda2, lambda3]
%            lambda(1): overall shrinkage for own lags
%            lambda(2): shrinkage for other variables' lags
%            lambda(3): shrinkage for intercept (constant)
%   sig2   - vector of error variances (length n) from OLS residuals or prior guess
%
% OUTPUT:
%   prior  - struct with fields:
%       beta0: prior means for VAR coefficients (size (n*(n*p+1)/n, n))
%       alp0: prior means for alphas (off-diagonal structural params)
%       Vbeta: prior variances for VAR coefficients
%       Valp: prior variances for alphas
%       nu: prior degrees of freedom for sigmas
%       S: prior scale parameters for sigmas

k_beta = n*(n*p+1);
k_alp = n*(n-1)/2;

% Initialize
prior.beta0 = zeros(k_beta/n, n);
prior.alp0 = zeros(k_alp, 1);
prior.Vbeta = zeros(k_beta/n, n);
prior.Valp = ones(k_alp,1) * 1e6; % diffuse prior variance for alphas (off-diagonal)
prior.nu = ones(n,1); % diffuse prior degrees of freedom
prior.S = sig2/2;     % prior scale set to half the variance (can adjust)

count_alp = 0;

for i = 1:n
    % Prior mean for coefficients: zeros except possibly intercept
    prior.beta0(:, i) = zeros(n*p+1,1);
    
    % Prior variance vector for beta coefficients
    Vbeta_vec = zeros(n*p+1,1);
    
    for j = 1:(n*p+1)
        if j == 1
            % Intercept prior variance
            Vbeta_vec(j) = lambda(3)^2;
        else
            lag_num = ceil((j-1)/n);        % lag number
            var_idx = mod(j-2, n) + 1;     % variable index
            
            if var_idx == i
                % Own lag
                Vbeta_vec(j) = (lambda(1)^2) / (lag_num^2 * sig2(var_idx));
            else
                % Other variables' lags
                Vbeta_vec(j) = (lambda(2)^2) / (lag_num^2 * sig2(var_idx));
            end
        end
    end
    
    prior.Vbeta(:, i) = Vbeta_vec;
    
    % Alphas prior mean and variance
    % Set prior mean zero and large variance for off-diagonal structural params
    prior.alp0(count_alp+1:count_alp+i-1) = zeros(i-1,1);
    prior.Valp(count_alp+1:count_alp+i-1) = ones(i-1,1) * 1e6;
    
    % Degrees of freedom and scale for error variance (diffuse)
    prior.nu(i) = 1;
    prior.S(i) = sig2(i)/2;
    
    count_alp = count_alp + i - 1;
end
