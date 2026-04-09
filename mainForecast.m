% This is the main run file for estimating the empirical forecasting application in
% Phella et al (2025)


clear all; clc;

addpath('Functions')

nloop = 20000;
burnin = 10000;

p = 4;   
data = xlsread('Final Dataset Forecast SSP.xlsx','hard_constraintapp1'); % load US data
[soft dat] = xlsread('Final Dataset Forecast SSP.xlsx','Soft_constraint I'); % load US data

hind = 1; % 0 - baseline & 1 - adverse
idx_ns = find(data(1,:)'==1);

data= data(2:end,:);
%data = data(100:end,:);

Y0 = data(1:p,:);  % save the first 8 obs as the initial conditions
Y = data(p+1:end,:);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
 [ml_opt,kappa] = get_OptKappa(Y0,Y,Z,p,[.04,.0016],'redu',idx_ns); 
sig2 = get_resid_var(Y0,Y,p);
prior_stru = prior_ACP_stru(n,p,kappa,sig2,idx_ns);
prior_redu = prior_ACP_redu(n,p,kappa,sig2,idx_ns);

% lambda=[0.1, 0.1, 1];
% prior_struc = prior_Minnesota_struc(n, p, lambda, sig2);
% prior_redu = prior_Minnesota_redu(n, p, lambda, sig2);
 
% tau   = 0.1;
% tailalpha = 0.5;
% prior_struc = prior_NG_struc(n,p,tau,tailalpha,sig2);
% prior_redu = prior_NG_redu(n,p,tau,tailalpha,sig2);

% prior_struc = prior_DL_struc(n,p,tau,sig2);
% prior_redu = prior_DL_redu(n,p,tau,sig2);

% prior_struc = prior_HS_struc(n, p, tau, sig2);
% prior_redu  = prior_HS_redu(n, p, tau, sig2);
%% process data
lag=p;T=T-lag;
%% prior
dim2=n*lag+1;dim=n*dim2; 

%% Set the conditions
T2=27; % Forecast horizon
n_con=1;% number of conditioning hard variables 
n_soft=2;%number of conditioning soft variables
n_t = n_con + n_soft;
M1=sparse(T2*n,T2*(n-n_con)); M2=sparse(T2*n,T2*n_con);
M1(find([ones(n-n_con,T2);zeros(n_con,T2)]==1),:)=speye(T2*(n-n_con));
M2(find([ones(n-n_con,T2);zeros(n_con,T2)]==0),:)=speye(T2*n_con);

M3=sparse(T2*(n-n_con),T2*(n-n_con-n_soft));M4=sparse(T2*(n-n_con),T2*n_soft);
M3(find([ones(n-n_con-n_soft,T2);zeros(n_soft,T2)]==1),:)=speye(T2*(n-n_con-n_soft));
M4(find([ones(n-n_con-n_soft,T2);zeros(n_soft,T2)]==0),:)=speye(T2*n_soft);
y_last=reshape(flip(Y(end-lag+1:end,:))',1,n*lag);
% if hind == 0
%     r_lower = soft(2:end,1)';
%     r_upper = soft(2:end,2)';
%     y_c = soft(2:end,3:4)';
% else
%     r_lower = soft(2:end,5)';
%     r_upper = soft(2:end,6)';
%    y_c = soft(2:end,7:8)';
% end

if hind == 0
    r_lower = soft(2:end,[1,3])';
    r_upper = soft(2:end,[2,4])';
    y_c = soft(2:end,5)';
else
    r_lower = soft(2:end,[6,8])';
    r_upper = soft(2:end,[7,9])';
    y_c = soft(2:end,10)';
end


Ind=cell(lag,1); 
for i=1:lag
    Ind{i,1}=spdiags(ones(T2-i,1),-i,T2,T2);
end

%% MCMC
Y_f=zeros(nloop-burnin,T2*(n-n_con));
Y_func=zeros(nloop-burnin,T2*n);
Y_fdiff=zeros(nloop-burnin,T2*(n-n_con));
Y_fdiffunc=zeros(nloop-burnin,T2*n);
Sens=zeros(T2*(n-n_t),T2,nloop-burnin);
invSig=speye(n);
Y_u=zeros(T2*(n-n_con),1);
store_invSig = zeros(n,n,nloop-burnin);
store_beta = zeros(dim2,n,nloop-burnin);

Y_new=Y;
Y0new = Y_new(1:p,:);
Ynew =  Y_new(p+1:end,:);

% MCMC starts here 
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC.... ');
disp(' ' );
start_time = clock;

for i=1:nloop

    %% Draw VAR and Sigma
    [alp,betahat,Sig] = sample_ThetaSig(Y0new,Ynew,p,prior_redu,1);
    [beta,Sighat] = getReducedForm(alp,betahat,Sig); 
    invSig = Sighat\speye(n);

    if i>burnin
    store_invSig(:,:,i-burnin) =  invSig;
      store_beta(:,:,i-burnin) = beta;



    end
  if ( mod( i, 2000 ) ==0 )
        disp(  [ num2str( i ) ' loops... ' ] )
    end 

end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

% Compute marginal likelihood
logML = logML_Chib_VAR(Y0new,Ynew,p,prior_redu,store_beta,store_invSig);

fprintf('Marginal likelihood (Chib) = %.3f\n', logML);

 %% Draw conditional forecast
for sim = 1:(nloop-burnin)

     beta = store_beta(:,:,sim)  ;
     invSig = store_invSig(:,:,sim);


     [alpha,B,K]=create_ABK2(y_last,beta,invSig,Ind,n,T2,lag);
     mu=alpha\B;

     Y_unconditional=mu+chol(K)\randn(T2*n,1);
     K_f=M1'*K*M1;      
     mu_u=K_f\(M1'*K*(mu-M2*y_c(:))); 
     V=M4'*mu_u;
      S=M4'/K_f*M4;
      Y_soft=mvrandn2(r_lower(:)-V,r_upper(:)-V,S,1)+V;
      K_hard=M3'*K_f*M3;      
      Y_hard=K_hard\(M3'*K_f*(mu_u-M4*Y_soft))+chol(K_hard)\randn(T2*(n-n_con-n_soft),1);
      Y_u=M4*Y_soft+M3*Y_hard;

      Y_f(sim,:)= Y_u;
      Y_func(sim,:) = Y_unconditional';
     Ycf =  reshape(Y_u,n-n_con,T2);
     Ynf =  reshape(Y_unconditional,n,T2);

      Y_fdiff(sim,:) = reshape(Ycf - Ynf(1:(n-n_con),:),T2*(n-n_con),1);

end

%%  conditional forecast  
base_forecast=reshape(mean(Y_f)',n-n_con,T2)'; % (n-1) x T*2
base_forecast16=reshape(quantile(Y_f,.16,1)',n-n_con,T2)'; % (n-1) x T*2
base_forecast84=reshape(quantile(Y_f,.84,1)',n-n_con,T2)'; % (n-1) x T*2

conforper = reshape(mean(Y_fdiff),n-n_con,T2)';
conforper16 = reshape(quantile(Y_fdiff,.16,1),n-n_con,T2)';
conforper84 = reshape(quantile(Y_fdiff,.84,1),n-n_con,T2)';

%%  unconditional forecast 
base_forecastunc=reshape(mean(Y_func)',n,T2)'; % (n-1) x T*2
base_forecastunc16=reshape(quantile(Y_func,.16,1)',n,T2)'; % (n-1) x T*2
base_forecastunc84=reshape(quantile(Y_func,.84,1)',n,T2)'; % (n-1) x T*2
conforperunc = reshape(mean(Y_fdiffunc)',T2,n);


lengthfor= length(base_forecastunc)-1;
%%

figure
subplot(2,2,[1,2])
plotx2r([repmat(data(end-lengthfor:end,2),1,3);[base_forecast(:,2) base_forecast16(:,2)  base_forecast84(:,2)] ],1:54)
hold on
plotx2n([repmat(data(end-lengthfor:end,2),1,3);[base_forecastunc(:,2) base_forecastunc16(:,2)  base_forecastunc84(:,2)]  ],1:54)
hold on
plot(1:T2,data(end-lengthfor:end,2),'k','LineWidth',3)
xticks([1 4 9 14 19 24 29 34 39 44 49 54])
xticklabels({'1997','2000','2005','2010','2015','2020','2025','2030','2035','2040','2045','2050'})
xlim([1 54])
grid on
ylabel('Levels')
legend('','Conditonal Forecast (Posterior Mean)','','Unconditional Forecast (Posterior Mean)','Observed Data','Actual Realization')
title('Temperatures','FontSize',20)
subplot(2,2,[3,4])
plotx2r([repmat(data(end-lengthfor:end,3),1,3);[base_forecast(:,3) base_forecast16(:,3)  base_forecast84(:,3)] ],1:54)
hold on
plotx2n([repmat(data(end-lengthfor:end,3),1,3);[base_forecastunc(:,3) base_forecastunc16(:,3)  base_forecastunc84(:,3)]  ],1:54)
hold on
plot(1:T2,data(end-lengthfor:end,3),'k','LineWidth',3)
xticks([1 4 9 14 19 24 29 34 39 44 49 54])
xticklabels({'1997','2000','2005','2010','2015','2020','2025','2030','2035','2040','2045','2050'})
xlim([1 54])
grid on
ylabel('Levels')
legend('','Conditonal Forecast (Posterior Mean)','','Unconditional Forecast (Posterior Mean)','Observed Data','Actual Realization')
title('Well-mixed GHG','FontSize',20)


% difference
figure
subplot(2,2,[1,2])
plotx2([zeros(1,3);[conforper(:,2) conforper16(:,2)  conforper84(:,2)  ]],1:T2+1)
grid on
xticks([1 8 13 18 23 28])
xticklabels({'2023','2030','2035','2040','2045','2050'})
xlim([1 28])
ylabel('Difference')
legend('','Conditional - Unconditional Forecast','','Actual Realization - Unconditional Forecast'),
title('Temperatures','FontSize',20)
subplot(2,2,[3,4])
plotx2([zeros(1,3);[conforper(:,3) conforper16(:,3)  conforper84(:,3)  ]],1:T2+1)
grid on
xticks([1 8 13 18 23 28])
xticklabels({'2023','2030','2035','2040','2045','2050'})
xlim([1 28])
ylabel('Difference')
legend('','Conditional - Unconditional Forecast','','Actual Realization - Unconditional Forecast'),
title('Well-mixed GHG','FontSize',20)






