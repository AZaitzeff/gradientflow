%nt=1024*2^6;
maxiter=100;
N=256;

nt=64;
tol=1e-4;
T=3/32;
eps=1/16;
[U,h]=initializecircle(N,eps);
Uinit=U;
error1=zeros(2,5);
dt=T/nt;
tic;
for t=1:nt
    curU=U;
    Un=U;
    for i=1:maxiter
        lastU=curU;
        b=constructb(curU,Un,N,dt,eps,h,1);
        func=@(x) applyoperator(x,curU,N,dt,eps,h,1);
        [x,flag]=pcg(func,b,1e-3,100,[],[],curU(:));
        %[x,flag]=minres(func,b,1e-6,100,[],[],curU(:));
        curU=curU+reshape(x, [N,N]);
        if (sum(abs(curU(:)-lastU(:)))/N)<tol
            break
        end
    end
    U=curU;
    
end
toc;