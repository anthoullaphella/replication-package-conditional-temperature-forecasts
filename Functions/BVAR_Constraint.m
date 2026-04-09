function [base_forecast, conforper, base_forecastunc, conforperunc] = BVAR_Constraint(Y0,Y,Z,r_lower, r_upper,y_c,n,p,T2,idx_ns)


lag=p;
nloop = 20000;
burnin = 10000;
dim2=n*lag+1;dim=n*dim2; 


Ind=cell(lag,1); 
for i=1:lag
    Ind{i,1}=spdiags(ones(T2-i,1),-i,T2,T2);
end

[ml_opt,kappa] = get_OptKappa(Y0,Y,Z,p,[.04,.0016],'redu',idx_ns); 
sig2 = get_resid_var(Y0,Y,p);
%prior_stru = prior_ACP_stru(n,p,kappa,sig2,idx_ns);
%prior_redu = prior_ACP_redu(n,p,kappa,sig2,idx_ns);

lambda=[0.1, 0.1, 1];
prior_struc = prior_Minnesota_struc(n, p, lambda, sig2);
prior_redu = prior_Minnesota_redu(n, p, lambda, sig2);


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

end