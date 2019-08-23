
N=256;

T=20;
eps=1/2;
[U,X,Y,h]=initializebigcircle(N,eps);
Uinit=U;

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
    Up=U;
    Lp=laplacian29(Up,N,h);
tic;
for t=1:nt
    if t==1
        L=laplacian29(U,N,h);
        U=U+dt*L-dt*1/eps^2*(2*U-6*U.^2+4*U.^3);
    else
        L=laplacian29(U,N,h);
        Utemp=U+3/2*(dt*L-dt*1/eps^2*(2*U-6*U.^2+4*U.^3))-1/2*(dt*Lp-dt*1/eps^2*(2*Up-6*Up.^2+4*Up.^3));
        Lp=L;
        Up=U;
        U=Utemp;
        
    end
    %F=U-dt*1/eps^2*(2*U-6*U.^2+4*U.^3
    
end
toc;
tU=trueU(1:8:end,1:8:end);
sqrt(h^2*sum((U(:)-tU(:)).^2))
%save(['results/feac' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
%error1(1,numt)=sqrt(h^2*sum((U(:)-trueU(:)).^2));
end
