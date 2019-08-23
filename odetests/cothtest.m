addpath('findgamma/')
[gamma,m]=getgamma(order,0);
T=2;
trueval=-2*acoth(exp(T)*coth(1));
Ns=[16,32,64,128,256,512];
maxiter=1000;

error1=zeros(2,6);
for iter=1:6
    N=Ns(iter);
    k=T/N;
    t=0:k:T;
    u=zeros(N,m+1);
    
    for i=1:N
        if i==1
          u(1,1)=-2;  
        else
          u(i,1)=u(i-1,m+1);
        end
        %u(i+1)=u(i)/(1+k);
        for j=2:m+1
            ubar=sum(gamma(j-1,1:j-1).*u(i,1:j-1));
            S=sum(gamma(j-1,1:j-1));
            curu=u(i,j-1);
            for z=1:maxiter
                curu=curu-(k*sinh(curu)+S*curu-ubar)/(k*cosh(curu)+S);
                if abs(-k*sinh(curu)+S*curu-ubar)<1e-14
                    break
                end
            end
            u(i,j)=curu;
        end
    end

    error1(1,iter)=abs(u(N,m+1)-trueval);
end
for i=2:6
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end