function y=applyoperator(x,ul,N,dt,eps,h,S)
    U=reshape(x, [N,N]);
    L=laplacian29(U,N,h);
    y=S*U-dt*L+dt/eps^2*(2-12*ul+12*ul.^2).*U;
    y=y(:);
end