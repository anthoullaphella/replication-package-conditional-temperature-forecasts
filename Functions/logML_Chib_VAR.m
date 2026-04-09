function logML = logML_Chib_VAR(Y0,Y,p,prior,store_beta,store_invSig)
% Computes the log marginal likelihood of a VAR using Chib (1995) method
%
% INPUTS:
%   Y0          - initial observations for lag construction (p x n)
%   Y           - estimation sample (T x n)
%   p           - number of lags
%   prior       - struct with prior parameters (beta0, Vbeta, alp0, Valp, nu, S)
%   store_beta  - posterior draws of VAR coefficients (dim x n x nsims)
%   store_invSig- posterior draws of inverse Sigma (n x n x nsims)
%
% OUTPUT:
%   logML       - log marginal likelihood of the VAR model

[T,n] = size(Y);
[dim2,~,nsims] = size(store_beta);
dim = dim2*n;

% Choose posterior evaluation point (posterior mean)
beta_star   = mean(store_beta,3);
invSig_star = mean(store_invSig,3);
Sig_star    = inv(invSig_star);

%% Log posterior at beta_star
% Prior contribution
Vbeta_inv = diag(prior.Vbeta(:).^(-1));
beta_vec = beta_star(:);
beta0_vec = prior.beta0(:);
log_prior = -0.5*( (beta_vec - beta0_vec)'*Vbeta_inv*(beta_vec - beta0_vec) );

% Likelihood contribution
log_like = 0;
% Construct lagged Y matrix including initial Y0
Y_full = [Y0(end-p+1:end,:); Y];
for t = p+1:T
    y_lag = reshape(flip(Y_full(t-p:t-1,:))', n*p, 1); % stacked lag vector
    X_t = [1; y_lag]; % intercept + lags
    mu_t = beta_star' * X_t;
    res_t = Y(t,:)' - mu_t;
    log_like = log_like - 0.5*( n*log(2*pi) + log(det(Sig_star)) + res_t'*(Sig_star\res_t) );
end

log_post = log_like + log_prior;

%% Estimate log posterior ordinate at beta_star
% Evaluate using Chib's method (average over draws)
log_post_draws = zeros(nsims,1);
for s = 1:nsims
    beta_s  = store_beta(:,:,s);
    invSig_s = store_invSig(:,:,s);
    Sig_s    = inv(invSig_s);

    log_like_s = 0;
    for t = p+1:T
        y_lag = reshape(flip(Y_full(t-p:t-1,:))', n*p, 1);
        X_t = [1; y_lag];
        mu_t = beta_s' * X_t;
        res_t = Y(t,:)' - mu_t;
        log_like_s = log_like_s - 0.5*( n*log(2*pi) + log(det(Sig_s)) + res_t'*(Sig_s\res_t) );
    end

    beta_vec_s = beta_s(:);
    log_prior_s = -0.5*( (beta_vec_s - beta0_vec)'*Vbeta_inv*(beta_vec_s - beta0_vec) );
    log_post_draws(s) = log_like_s + log_prior_s;
end

% log posterior ordinate at beta_star using log-sum-exp
m = max(log_post_draws);
log_post_ord = m + log(mean(exp(log_post_draws - m)));

%% Log marginal likelihood
logML = log_post - log_post_ord;

end