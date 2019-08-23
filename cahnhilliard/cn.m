addpath('../findgamma/')
maxiter=200;
% fac=1024/N;
% 
% vars=load('feulerac');
% trueU=vars.U(1:fac:end,1:fac:end);


tol=1e-11;
tol1=1e-2;
T=20;
eps=1/2;

nts=[2^10];


for numt=1:1
    N=256;
    [Uinit,X,Y,h]=initializebigoval(N,eps);
    force=1/(4*pi)*exp(-(X.^2+Y.^2)/4);
    nt=nts(numt);
    U=Uinit;
    dt=T/nt;
    tic;
for t=1:nt
    curU=U;
    Un=U;
    
        
    for i=1:maxiter
        lastU=curU;
        b=constructbforcn(curU,Un,N,dt,eps,h);
        func=@(x) 1/2*applyoperator(x,curU,N,dt,eps,h,2);
        [x,flag,res,iter]=pcg(func,b,tol1,400);
        %[x,flag]=minres(func,b,1e-6,100,[],[],curU(:));
        curU=curU+reshape(x, [N,N]);
        if (max(abs(curU(:)-lastU(:))))<tol
            break
        end
        
        %val=normdiff(curU,Un,N,dt,eps,h,S,force)
    end
    U=curU;
    
end
toc;
%save(['results/cnac' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
%sqrt(sum((U(:)-trueU(:)).^2))/sqrt(sum((Uinit(:)-trueU(:)).^2))
%error1(1,numt)=sqrt(sum((U(:)-trueU(:)).^2))/sqrt(sum(trueU(:).^2));
end