function prior = prior_NG_redu(n,p,tau,alpha,sig2)
% Normal–Gamma prior on reduced-form VAR parameters

k_beta = n*(n*p+1);

prior_stru = prior_NG_struc(n,p,tau,alpha,sig2);

prior.beta0 = prior_stru.beta0;
prior.alp0  = prior_stru.alp0;
prior.Valp  = prior_stru.Valp;
prior.nu    = prior_stru.nu;
prior.S     = prior_stru.S;

prior.Vbeta = zeros(k_beta/n,n);

for ii = 1:n
    for jj = 1:(n*p+1)
        if ii == 1
            prior.Vbeta(jj,ii) = prior_stru.Vbeta(jj,ii);
        else
            prior.Vbeta(jj,ii) = prior_stru.Vbeta(jj,ii) + ...
                sum(prior_stru.Vbeta(jj,1:ii-1) + ...
                prior_stru.beta0(jj,1:ii-1).^2 ./ sig2(1:ii-1)');
        end
    end
end

end