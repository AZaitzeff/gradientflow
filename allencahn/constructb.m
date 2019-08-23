function b=constructb(u,un,N,dt,eps,h,S,force)
    L=laplacian29(u,N,h);
    b=S*u-un-dt*L+dt/eps^2*(2*u-6*u.^2+4*u.^3)-dt*(force);
    b=-b(:);
end