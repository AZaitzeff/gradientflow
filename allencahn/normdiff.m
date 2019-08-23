function val=normdiff(u,un,N,dt,eps,h,S,force)
    L=laplacian29(u,N,h);
    
    y=S*u-dt*L+dt/eps^2*(2*u-6*u.^2+4*u.^3)-dt*(force-u);
    %L=(4/3*(U(2:N+1,3:N+2)+U(4:N+3,3:N+2)+U(3:N+2,2:N+1)+U(3:N+2,4:N+3))-...
    %    1/12*(U(1:N,3:N+2)+U(5:N+4,3:N+2)+U(3:N+2,1:N)+U(3:N+2,5:N+4)))/h^2;
    
    %y=U(3:N+2,3:N+2)-dt*L;
    val=mean(abs(y(:)-un(:)));
end