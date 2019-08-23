function U=semiimplicitsteps(U,invmatrix,N,h,dt,eps,nt)

for t=1:nt
    L=laplacian9(U,N,h);
    F=U+dt/2*L-dt*1/eps^2*(2*U-6*U.^2+4*U.^3)-...
        dt^2/2*(1/eps^2*(2-12*U+12*U.^2).*(L-1/eps^2*(2*U-6*U.^2+4*U.^3)));
    %F=U-dt*1/eps^2*(2*U-6*U.^2+4*U.^3);
    Fbar=fft2(F);
    Ubar=Fbar./(invmatrix);
    U=real(ifft2(Ubar));
    
end