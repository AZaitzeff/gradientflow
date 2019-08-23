function L=laplacian5(U,N,h)


    U=padarray(U,[2,2],'circular');
    L=(-4*U(3:N+2,3:N+2)+1*(U(2:N+1,3:N+2)+U(4:N+3,3:N+2)+U(3:N+2,2:N+1)+U(3:N+2,4:N+3)))/h^2;