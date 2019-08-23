nt=2^10;
N=512;
T=1/2;
[Uinit,trueU,X,Y,h]=initializedata2d(N,T);


l=(0:(N-1));
k=(0:(N-1));

%invmatrix=1-dt/(2*h^2)*(-4976/1152+3/2*(cos(2*pi*l/N)+cos(2*pi*k/N))...
%    +3/40*(cos(4*pi*l/N)+cos(4*pi*k/N))+1/45*(cos(6*pi*l/N)+cos(6*pi*k/N))...
%    +2*(cos(2*pi*l/N)*cos(2*pi*k/N))+1/8*(cos(4*pi*l/N)*cos(4*pi*k/N))...
%    -1/2*(cos(4*pi*l/N)*cos(2*pi*k/N)+cos(2*pi*l/N)*cos(4*pi*k/N)));
nts=[2^6,2^7,2^8,2^9,2^10,2^11,2^12,2^13,2^14];
NT=numel(nts);
error1=zeros(2,NT);
error2=zeros(2,NT);
error3=zeros(2,NT);
for numt=1:1
    U=Uinit;
    nt=nts(numt);
    dt=T/nt;
    invmatrix=1-dt/(2*h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
    tic;
for t=1:nt
    L=laplacian9(U,N,h);
    F=U+dt/2*L;
    %F=U-dt*1/eps^2*(2*U-6*U.^2+4*U.^3);
    Fbar=fft(F);
    Ubar=Fbar./(invmatrix);
    U=real(ifft(Ubar));
end
toc;
error1(1,numt)=sqrt(h*sum((U(:)-trueU(:)).^2));
error2(1,numt)=sqrt(h*sum((U(:)-trueU(:)).^2))/sqrt(h*sum((trueU(:)).^2));
error3(1,numt)=max(abs(U(:)-trueU(:)));
end


for i=2:NT
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error2(2,i)=(log(error2(1,i-1))-log(error2(1,i)))/log(2);
    error3(2,i)=(log(error3(1,i-1))-log(error3(1,i)))/log(2);
end