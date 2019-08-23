function D=builddlap(N,h)
    D=sparse(N^2,N^2);
    vec=[[1,2,-1,-2,0,0,0,0];[0,0,0,0,1,2,-1,-2]];
    fac=[4/3,-1/12,4/3,-1/12,4/3,-1/12,4/3,-1/12];
    for n=1:N^2
        D(n,n)=-5/h^2;
        [i,j] = ind2sub([N N],n);
        for k=1:8
            ni=mod(i-1+vec(1,k),N)+1;
            nj=mod(j-1+vec(1,k),N)+1;
            nn=sub2ind([N N],[ni nj]);
            D(n,nn)=fac(k)/h^2;
        end
        
        
    end



end