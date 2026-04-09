function logp = logmvnpdf(y,mu,Sigma)
% LOGMVNPDF  Log-density of multivariate normal N(mu,Sigma)
%
% INPUTS:
%   y     : n x 1 vector
%   mu    : n x 1 mean vector
%   Sigma : n x n covariance matrix
%
% OUTPUT:
%   logp  : scalar log-density

n = length(y);

% Cholesky factorization (Sigma must be SPD)
C = chol(Sigma,'lower');

% Solve for standardized residual
z = C \ (y - mu);

% Log determinant via Cholesky
logdetSigma = 2 * sum(log(diag(C)));

% Log-density
logp = -0.5 * ( n*log(2*pi) + logdetSigma + z'*z );

end