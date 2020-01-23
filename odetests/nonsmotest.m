addpath('../findgamma/')
order=3;
[gamma,m]=getgamma(order,0);
%gamma=[[1]];
%m=1;
T=2.5;
trueval=-1/2*(T-1)+1;
Ns=2.^(4:21);
n=numel(Ns);
maxiter=1000;
dts=zeros(1,n);
error1=zeros(1,n);
for iter=1:1
    N=Ns(iter);
    k=T/N;
    dts(iter)=k;
    t=0:k:T;
    u=zeros(N,m+1);
    
    for i=1:N
        if i==1
          u(1,1)=2;  
        else
          u(i,1)=u(i-1,m+1);
        end
        %u(i+1)=u(i)/(1+k);
        for j=2:m+1
            ubar=sum(gamma(j-1,1:j-1).*u(i,1:j-1));
            S=sum(gamma(j-1,1:j-1));
            curu=solveimplpwcons(k/S,ubar/S);
            %curu=solveimplpwcons2(S/k,ubar/S);
            u(i,j)=curu;
        end
    end

    error1(iter)=abs(u(N,m+1)-trueval);
end
% for i=2:6
%     error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
% end