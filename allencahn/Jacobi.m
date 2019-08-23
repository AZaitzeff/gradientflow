function U=Jacobi(U,ul,N,dt,eps,h,S,b,nt)
    U=reshape(U,[N,N]);
    b=reshape(b,[N,N]);
    for i=1:nt
        U=padarray(U,[2,2],'circular');
        R=-dt*(4/3*(U(2:N+1,3:N+2)+U(4:N+3,3:N+2)+U(3:N+2,2:N+1)+U(3:N+2,4:N+3))-...
            1/12*(U(1:N,3:N+2)+U(5:N+4,3:N+2)+U(3:N+2,1:N)+U(3:N+2,5:N+4)))/h^2;
        D=S+dt*(5)/h^2+dt/eps^2*(2-12*ul+12*ul.^2);
        U=(b-R)./D;
    end
end