function L=bilaplacian71d(U,N,h)

    U=padarray(U,[0,3],'circular');
    L=(28/3*U(1,4:N+3)...
        -13/2*(U(1,3:N+2)+U(1,5:N+4))...
        +2*(U(1,2:N+1)+U(1,6:N+5))...
        -1/6*(U(1,1:N)+U(1,7:N+6)))/h^4;

