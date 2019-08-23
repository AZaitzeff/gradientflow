function y=applyoperator(x,ul,N,dt,eps,h,S)
    U=reshape(x, [N,N]);
    L=laplacian29(U,N,h);
    LW=laplacian29(L-1/eps^2*(2-12*ul+12*ul.^2).*U,N,h);
    y=S*U+dt*LW;
    y=y(:);
end