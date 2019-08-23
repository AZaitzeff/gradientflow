function L=laplacian9s(U,N,h)


    U=padarray(U,[2,2],'circular');
    L=(-10/3*U(3:N+2,3:N+2)+2/3*(U(2:N+1,3:N+2)+U(4:N+3,3:N+2)+U(3:N+2,2:N+1)+U(3:N+2,4:N+3))+...
        1/6*(U(2:N+1,2:N+1)+U(4:N+3,4:N+3)+U(4:N+3,2:N+1)+U(2:N+1,4:N+3)))/h^2;