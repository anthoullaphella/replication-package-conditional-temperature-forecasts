function prior = prior_HS_redu(n, p, tau, sig2)
% Constructs a Horseshoe (global shrinkage) prior on reduced-form VAR parameters.
%
% This function first builds the Horseshoe prior on the structural
% parameterization and then maps it into the reduced-form VAR
% parameterization following Chan (2021).
%
% INPUTS:
%   n     - number of variables
%   p     - number of lags
%   tau   - global shrinkage parameter (scalar)
%   sig2  - vector of error variances (length n)
%
% OUTPUT:
%   prior - struct with fields:
%       beta0  : prior means for VAR coefficients
%       alp0   : prior means for structural parameters
%       Vbeta  : prior variances for VAR coefficients
%       Valp   : prior variances for structural parameters
%       nu     : prior degrees of freedom for sigmas
%       S      : prior scale parameters for sigmas



% -------------------------------------------------------------------------
% Build structural Horseshoe prior
% -------------------------------------------------------------------------

prior_stru = prior_HS_struc(n, p, tau, sig2);

% Number of coefficients per equation
k_beta = n * p + 1;

% -------------------------------------------------------------------------
% Copy over prior means and variance parameters
% -------------------------------------------------------------------------

prior.beta0 = prior_stru.beta0;
prior.alp0  = prior_stru.alp0;
prior.Valp  = prior_stru.Valp;

prior.nu = prior_stru.nu;
prior.S  = prior_stru.S;

% -------------------------------------------------------------------------
% Map structural variances to reduced-form variances
% -------------------------------------------------------------------------

prior.Vbeta = zeros(k_beta, n);

for ii = 1:n
    for jj = 1:k_beta

        if ii == 1
            % First equation unchanged
            prior.Vbeta(jj,ii) = prior_stru.Vbeta(jj,ii);
        else
            % Accumulate structural uncertainty from previous equations
            prior.Vbeta(jj,ii) = prior_stru.Vbeta(jj,ii) ...
                + sum( prior_stru.Vbeta(jj,1:ii-1) ...
                + prior_stru.beta0(jj,1:ii-1).^2 ./ sig2(1:ii-1)' );
        end

    end
end

end