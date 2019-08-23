addpath('findgamma/')
[gamma,m]=getgamma(order,0);

T=2;
Ns=[16,32,64,128,256,512];

error1=zeros(2,6);
for iter=1:6
    N=Ns(iter);
    k=T/N;
    t=0:k:T;
    u=zeros(N,m+1);
    func=@(x) x^2/2;
    funcx=@(x) x;
    
    for i=1:N
        if i==1
          u(1,1)=1;  
        else
          u(i,1)=u(i-1,m+1);
        end
        %u(i+1)=u(i)/(1+k);
        for j=2:m+1
            u(i,j)=sum(gamma(j-1,1:j-1).*u(i,1:j-1))/(sum(gamma(j-1,1:j-1))+k);
        end
    end

    error1(1,iter)=abs(u(N,m+1)-exp(-T));
end

for i=2:6
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end