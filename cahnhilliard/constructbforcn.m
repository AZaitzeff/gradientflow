function b=constructbforcn(u,un,N,dt,eps,h)
    L=laplacian9((u+un)/2,N,h);
    b=u-un+laplacian9(dt*L-1/2*dt/eps^2*(2*u-6*u.^2+4*u.^3)-1/2*dt/eps^2*(2*un-6*un.^2+4*un.^3),N,h);
    b=-b(:);
end