addpath('findgamma/')
[gamma,m]=getgamma(order,0);
maxiter=200;


tol=1e-9;
tol1=1e-3;
T=20;
eps=1;

nts=[4,2^3,2^4,2^5,2^6];


for numt=1:5
    N=ceil(512*2^((numt-1)/4));
    Uall=zeros(m+1,N,N);
    [Uinit,X,Y,h]=initializebigoval(N,eps);
    nt=nts(numt);
    U=Uinit;
    dt=T/nt;
    tic;
for t=1:nt
    curU=U;
    Uall(1,:,:)=U;
    for step=1:m
        Un=squeeze(Uall(1,:,:)*gamma(step,1));
        for zm=2:step
            Un=Un+squeeze(Uall(zm,:,:)*gamma(step,zm));
        end
        S=sum(gamma(step,:));
        
        for i=1:maxiter
            lastU=curU;
            b=constructb(curU,Un,N,dt,eps,h,S);
            func=@(x) applyoperator(x,curU,N,dt,eps,h,S);
            [x,flag,res,iter]=pcg(func,b,tol1,30000);
            %[x,flag]=minres(func,b,1e-6,100,[],[],curU(:));
            curU=curU+reshape(x, [N,N]);
            if (max(abs(curU(:)-lastU(:))))<tol
                break
            end
        end
        %val=normdiff(curU,Un,N,dt,eps,h,S,force)
        Uall(step+1,:,:)=curU;
    end
    U=curU;
    
end
toc;
save(['results/multistepch' num2str(order) 's' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
%sqrt(sum((U(:)-trueU(:)).^2))/sqrt(sum((Uinit(:)-trueU(:)).^2))
%error1(1,numt)=sqrt(sum((U(:)-trueU(:)).^2))/sqrt(sum(trueU(:).^2));
end