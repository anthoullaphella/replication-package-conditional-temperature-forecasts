function prior = prior_DL_struc(n,p,tau,sig2)
% Constructs Dirichlet–Laplace-style prior on structural VAR parameters
%
% INPUTS:
%   n    - number of variables
%   p    - number of lags
%   tau  - global shrinkage parameter (scalar, e.g. 0.1)
%   sig2 - reduced-form residual variances (n x 1)
%
% OUTPUT:
%   prior - struct compatible with sample_ThetaSig

k_beta = n*(n*p+1);
k_alp  = n*(n-1)/2;

prior.beta0 = zeros(k_beta/n,n);   % zero means
prior.alp0  = zeros(k_alp,1);

prior.Vbeta = zeros(k_beta/n,n);
prior.Valp  = zeros(k_alp,1);

prior.nu = ones(n,1);
prior.S  = sig2/2;

% ---- Dirichlet weights (fixed expectation) ----
K = n*p + 1;
a = 1/K;                     % sparsity level
phi = ones(K,1)/K;           % E[Dirichlet(a)]

count_alp = 0;

for ii = 1:n
    % VAR coefficients
    for j = 1:K
        prior.Vbeta(j,ii) = 2 * tau^2 * phi(j)^2 / sig2(ii);
    end

    % Structural coefficients (lower triangular)
    for j = 1:(ii-1)
        prior.Valp(count_alp+j) = 2 * tau^2 / sig2(j);
    end

    count_alp = count_alp + ii - 1;
end

end