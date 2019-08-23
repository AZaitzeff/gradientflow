nt=2^20;
N=4096;
T=20;
eps=1/2;
[U,X,Y,h]=initializebigcircle(N,eps);
Uinit=U;
dt=T/nt;
[l,k]=meshgrid(0:(N-1));
invmatrix=1-dt/(2*h^2)*(-5+8/3*cos(2*pi*l/N)+8/3*cos(2*pi*k/N)-1/6*cos(4*pi*l/N)-1/6*cos(4*pi*k/N));
%invmatrix=1-dt/(2*h^2)*(-4976/1152+3/2*(cos(2*pi*l/N)+cos(2*pi*k/N))...
%    +3/40*(cos(4*pi*l/N)+cos(4*pi*k/N))+1/45*(cos(6*pi*l/N)+cos(6*pi*k/N))...
%    +2*(cos(2*pi*l/N)*cos(2*pi*k/N))+1/8*(cos(4*pi*l/N)*cos(4*pi*k/N))...
%    -1/2*(cos(4*pi*l/N)*cos(2*pi*k/N)+cos(2*pi*l/N)*cos(4*pi*k/N)));
tic;
for t=1:nt

    L=laplacian9(U,N,h);
    F=U+dt/2*L-dt*1/eps^2*(2*U-6*U.^2+4*U.^3)-1/eps^2*dt^2/2*(2-12*U+12*U.^2).*(L-1/eps^2*(2*U-6*U.^2+4*U.^3));
    Fbar=fft2(F);
    Ubar=Fbar./(invmatrix);
    U=real(ifft2(Ubar));
    %F=U-dt*1/eps^2*(2*U-6*U.^2+4*U.^3);
    
    
end
toc;
%tU=trueU(1:2:end,1:2:end);
%sqrt(h^2*sum((U(:)-tU(:)).^2))
save(['results/semieulerac' num2str(nt) 'N' num2str(N)],'T','dt','N','U','Uinit','eps')