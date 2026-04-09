% This is the main run file for estimating the empirical counterfactual forecasting application in
% Phella et al (2025)
addpath('Functions')

clear all; clc;
nloop  = 20000;
burnin = 10000;

p = 4;   
T_start = p + 156;

data_full = xlsread('Final Dataset Forecast SSP Counterfactual.xlsx','hard_constraintapp1');
data_full = data_full(2:end,:);

[Tfull,n] = size(data_full);

%% ============ STORAGE ===================
logPL = zeros(Tfull - T_start,1);
cnt = 1;
%% ========================================

for t = T_start:(Tfull-1)

    fprintf('Forecast origin t = %d\n',t);

    %% ---------- Estimation sample ----------
    data = data_full(1:t,:);

    Y0 = data(1:p,:);
    Y  = data(p+1:end,:);
    [T,~] = size(Y);

    tmpY = [Y0(end-p+1:end,:); Y];
    Z = zeros(T,n*p); 
    for ii=1:p
        Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
    end
    Z = [ones(T,1) Z];

    %% ---------- Prior ----------
    sig2 = get_resid_var(Y0,Y,p);
    [~,kappa] = get_OptKappa(Y0,Y,Z,p,[.04,.0016],'redu',[]);
    prior_redu = prior_ACP_redu(n,p,kappa,sig2,[]);

    %% ---------- MCMC ----------
    store_beta = zeros(n*p+1,n,nloop-burnin);
    store_Sig  = zeros(n,n,nloop-burnin);

    for i = 1:nloop

        [alp,betahat,Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,1);
        [beta,Sighat] = getReducedForm(alp,betahat,Sig);

        if i > burnin
            store_beta(:,:,i-burnin) = beta;
            store_Sig(:,:,i-burnin)  = Sighat;
        end
    end

    %% ---------- Predictive likelihood ----------
    y_next = data_full(t+1,:)';

    y_lag = reshape(flip(data(end-p+1:end,:))',n*p,1);
    X1 = [1; y_lag];

    logdens = zeros(nloop-burnin,1);

    for s = 1:(nloop-burnin)
        beta = store_beta(:,:,s);
        Sig  = store_Sig(:,:,s);
        mu   = beta' * X1;

        logdens(s) = logmvnpdf(y_next,mu,Sig);
    end

    %% ---------- Log-sum-exp ----------
    m = max(logdens);
    logPL(cnt) = m + log(mean(exp(logdens - m)));
    cnt = cnt + 1;

end

%% ============ FINAL SCORE =================
CLPS = sum(logPL);
fprintf('\nCumulative log predictive score = %.3f\n',CLPS);