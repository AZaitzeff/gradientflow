%nt=1024*2^6;
vars=load(['results/feacpolar67108864N4096']);
%vars=load(['results/multistepac' num2str(2^10) 'N' num2str(1024)]);

eps=1/2;
%[~,r,~]=initializebigcirclepolar(4096,eps);
nts=[2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10];%2^9];
n=numel(nts);
error1=zeros(2,n);
error2=zeros(2,n);
error3=zeros(2,n);
for numt=1:4
    nt=nts(numt);
    %N=ceil(512*2^((numt-2)/2));
    N=256;
    [~,X,Y,h]=initializebigcircle(N,1);
    load(['results/cnac' num2str(nt) 'N' num2str(N)]);
    trueU=zeros(N,N);
    R=sqrt(X.^2+Y.^2);
    mask=(R)>14.9;
    trueU=vars.U;
    trueU(mask)=0;
    %trueU = interp2(tX,tY,vars.U,X,Y, 'spline');
    trueU(~mask) = interp1(r,vars.U,R(~mask), 'spline');
    error1(1,numt)=sqrt(h^2*sum((U(:)-trueU(:)).^2));
    error2(1,numt)=max(abs(U(:)-trueU(:)));
    error3(1,numt)=h^2*sum(abs(U(:)-trueU(:)));
end

for i=2:n
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error2(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error3(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end