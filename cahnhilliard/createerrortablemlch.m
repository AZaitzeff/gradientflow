
vars=load(['results/siech4194304N1024.mat']);
trueU=vars.U;
eps=1;
[~,trueX,trueY,h]=initializebigoval(1024,eps);
%trueU=vars.U(1:end,1:end);

nts=[4,2^3,2^4,2^5,2^6];
%nts=[2^8,2^9,2^10,2^11,2^12,2^13];
n=numel(nts);
error1=zeros(2,n);
error2=zeros(2,n);
error3=zeros(2,n);
for numt=1:n
    nt=nts(numt);
    N=ceil(512*2^((numt-1)/4));
    %N=128;
    [~,X,Y,h]=initializebigoval(N,eps);
    load(['results/multistepch' num2str(order) 's' num2str(nt) 'N' num2str(N)]);
    tU= interp2(trueX,trueY,trueU,X,Y,'spline');
    error1(1,numt)=sqrt(h^2*sum((U(:)-tU(:)).^2));
    error2(1,numt)=max(abs(U(:)-tU(:)));
    error3(1,numt)=h^2*sum(abs(U(:)-tU(:)));
end

for i=2:n
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error2(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error3(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end