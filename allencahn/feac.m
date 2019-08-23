
N=256;

T=20;
eps=1/2;
[U,X,Y,h]=initializebigcircle(N,eps);
Uinit=U;
dt=T/nt;
%force=(erf((5-X))+erf((5+X))).*(erf((5-Y))+erf((5+Y)))/4;
%invmatrix=1-dt/(2*h^2)*(-4976/1152+3/2*(cos(2*pi*l/N)+cos(2*pi*k/N))...
%    +3/40*(cos(4*pi*l/N)+cos(4*pi*k/N))+1/45*(cos(6*pi*l/N)+cos(6*pi*k/N))...
%    +2*(cos(2*pi*l/N)*cos(2*pi*k/N))+1/8*(cos(4*pi*l/N)*cos(4*pi*k/N))...
%    -1/2*(cos(4*pi*l/N)*cos(2*pi*k/N)+cos(2*pi*l/N)*cos(4*pi*k/N)));
nts=[2^16];

error1=zeros(2,4);
for numt=1:1
    nt=nts(numt);
dt=T/nt;

    U=Uinit;
tic;
for t=1:nt
    L=laplacian29(U,N,h);
    U=U+dt*L-dt*1/eps^2*(2*U-6*U.^2+4*U.^3);
    %F=U-dt*1/eps^2*(2*U-6*U.^2+4*U.^3
    
end
toc;
%save(['results/feac' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
%error1(1,numt)=sqrt(h^2*sum((U(:)-trueU(:)).^2));
end

%for i=2:4
%    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
%end
%save('feulerac','T','dt','N','U','eps')