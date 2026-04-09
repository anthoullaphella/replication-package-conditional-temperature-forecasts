function prior = prior_NG_struc(n,p,tau,alpha,sig2)
% Constructs Normal–Gamma prior on structural VAR parameters
%
% INPUTS:
%   n     - number of variables
%   p     - number of lags
%   tau   - global shrinkage parameter
%   alpha - shape parameter of Gamma mixing
%   sig2  - residual variances (n x 1)
%
% OUTPUT:
%   prior - struct compatible with sample_ThetaSig

k_beta = n*(n*p+1);
k_alp  = n*(n-1)/2;

prior.beta0 = zeros(k_beta/n,n);
prior.alp0  = zeros(k_alp,1);

prior.Vbeta = zeros(k_beta/n,n);
prior.Valp  = zeros(k_alp,1);

prior.nu = ones(n,1);
prior.S  = sig2/2;

% Expected local variance under Gamma
E_lambda = 1/alpha;

count_alp = 0;

for ii = 1:n
    % VAR coefficients
    for j = 1:(n*p+1)
        prior.Vbeta(j,ii) = tau^2 * E_lambda / sig2(ii);
    end

    % Structural coefficients
    for j = 1:(ii-1)
        prior.Valp(count_alp+j) = tau^2 * E_lambda / sig2(j);
    end

    count_alp = count_alp + ii - 1;
end

end