function b=constructb(u,un,N,dt,eps,h,S)
    L=laplacian29(u,N,h);
    b=S*u-un+dt*laplacian29(L-1/eps^2*(2*u-6*u.^2+4*u.^3),N,h);
    b=-b(:);
end