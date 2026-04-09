function prior = prior_HS_struc(n, p, tau, sig2)% Constructs a Horseshoe (global shrinkage) prior on structural VAR parameters.
%
% This is a Gaussian scale-mixture approximation to the Horseshoe prior,
% implemented in a conjugate form so it is compatible with
% sample_ThetaSig.m (no sampler changes required).
%
% The prior takes the form:
%   beta_j | sigma_i^2 ~ N(0, sigma_i^2 * tau^2 * lambda_j^2)
%
% where:
%   tau      - global shrinkage parameter
%   lambda_j - local (fixed) shrinkage weights by coefficient type
%
% INPUTS:
%   n     - number of variables
%   p     - number of lags
%   tau   - global shrinkage parameter (scalar, e.g. 0.1)
%   sig2  - vector of error variances (length n)
%
% OUTPUT:
%   prior - struct with fields:
%       beta0  : prior means for VAR coefficients (zeros)
%       alp0   : prior means for structural parameters (zeros)
%       Vbeta  : prior variances for VAR coefficients
%       Valp   : prior variances for structural parameters
%       nu     : prior degrees of freedom for sigmas
%       S      : prior scale parameters for sigmas



% Dimensions
k_beta = n*p + 1;
k_alp  = n*(n-1)/2;

% Zero-mean shrinkage
prior.beta0 = zeros(k_beta,n);
prior.alp0  = zeros(k_alp,1);

% Allocate variances
prior.Vbeta = zeros(k_beta,n);
prior.Valp  = zeros(k_alp,1);

% Inverse-Gamma prior for sigmas
prior.nu = 2*ones(n,1);
prior.S  = sig2(:);

count_alp = 0;

for i = 1:n
    for j = 1:k_beta

        % Identify lag and variable
        if j == 1
            lambda2 = 10;     % intercept
        else
            lag = ceil((j-1)/n);
            var = mod(j-1,n); if var==0, var=n; end

            if var==i && lag==1
                lambda2 = 5;      % own first lag
            elseif var==i
                lambda2 = 1;      % own higher lag
            else
                lambda2 = 0.1;    % cross-variable lag
            end
        end

        % Horseshoe-style variance
        prior.Vbeta(j,i) = tau^2 * lambda2 / sig2(i);
    end

    % Structural parameters: aggressive shrinkage
    for j = 1:(i-1)
        prior.Valp(count_alp+j) = tau^2 / sig2(i);
    end

    count_alp = count_alp + i - 1;
end

end