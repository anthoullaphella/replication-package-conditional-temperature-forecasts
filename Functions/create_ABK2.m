function [alpha,B,K]=create_ABK2(y_last,beta,invS,Ind,n,T2,lag)


 alpha=speye(n*T2);
 B=repmat(beta(1,:)',T2,1);
   K=kron(speye(T2),invS);
M=invS*beta(2:end,:)';
M2=beta(2:end,:)*M;

 for i=1:lag
     indx=(i-1)*n+1:i*n;
     alpha=alpha-kron(Ind{i,1},beta(indx+1,:)');
     B(indx)=B(indx)+(y_last(1:end-(i-1)*n)*beta((i-1)*n+2:end,:))';
      K=K-kron(Ind{i,1},M(:,indx))-kron(Ind{i,1}',M(:,indx)');

for j=1:lag
    indjx=(j-1)*n+1:j*n; 
    K=K+kron(Ind{i,1}'*Ind{j,1},M2(indx,indjx));
end
 end

        