
vars=load(['results/si2ac262144N2048.mat']);
trueU=vars.U;
eps=1;
nts=[2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10];
n=numel(nts);
error1=zeros(2,n);
error2=zeros(2,n);
error3=zeros(2,n);
for numt=1:n
    nt=nts(numt);
    %N=ceil(512*2^((numt-2)/2));
    N=2048;
    [~,X,Y,h]=initializebigcircle(N,1);
    load(['results/multistepac' num2str(order) 's' num2str(nt) 'N' num2str(N)]);
    error1(1,numt)=sqrt(h^2*sum((U(:)-trueU(:)).^2));
    error2(1,numt)=max(abs(U(:)-trueU(:)));
    error3(1,numt)=h^2*sum(abs(U(:)-trueU(:)));
end

for i=2:n
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error2(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error3(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end