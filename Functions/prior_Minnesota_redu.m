function prior = prior_Minnesota_redu(n, p, lambda, sig2)
% Constructs Minnesota prior on reduced-form VAR parameters and
% transforms it to implied structural prior, matching original structure.
%
% INPUTS:
%   n      - number of variables
%   p      - number of lags
%   lambda - vector of hyperparameters [lambda1, lambda2, lambda3]
%   sig2   - vector of residual variances (length n)
%
% OUTPUT:
%   prior - struct with fields:
%       beta0, alp0, Vbeta, Valp, nu, S
%
% This uses the Minnesota prior on reduced-form coefficients, then
% derives the implied prior on structural parameters.

% Get Minnesota prior on structural form first (to get alp0, Valp, nu, S)
prior_stru = prior_Minnesota_struc(n, p, lambda, sig2);

k_beta = n*(n*p+1);

% Initialize reduced-form prior
prior.beta0 = zeros(k_beta/n, n);
prior.Vbeta = zeros(k_beta/n, n);

prior.alp0 = prior_stru.alp0;
prior.Valp = prior_stru.Valp;
prior.nu = prior_stru.nu;
prior.S = prior_stru.S;

for ii = 1:n
    for jj = 1:(n*p+1)
        if ii == 1
            % For first variable, reduced-form prior variance = structural prior variance
            prior.Vbeta(jj, ii) = prior_stru.Vbeta(jj, ii);
            prior.beta0(jj, ii) = prior_stru.beta0(jj, ii);
        else
            % For variables 2,...,n, add structural prior variances and
            % contributions from earlier variables
            prior.Vbeta(jj, ii) = prior_stru.Vbeta(jj, ii) + ...
                sum(prior_stru.Vbeta(jj, 1:ii-1) + (prior_stru.beta0(jj, 1:ii-1).^2)./sig2(1:ii-1)');
            
            prior.beta0(jj, ii) = prior_stru.beta0(jj, ii); % remains zero here
        end
    end
end

end