addpath('findgamma/');

if order==2
[gamma,m]=getgamma(2,0);
maxiter=200;

tol=1e-9;
tol1=1e-3;
T=20;
eps=1;
nts=[2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10];

for numt=1:10
    %N=ceil(512*2^((numt-2)/2));
    N=2048;
    Uall=zeros(m+1,N,N);
    [Uinit,X,Y,h]=initializebigcircle(N,eps);
    %force=1/4*1./(1+exp(-1/eps^2*(4-sqrt(X.^2+Y.^2))));
    force=0;
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
            b=constructb(curU,Un,N,dt,eps,h,S,force);
            func=@(x) applyoperator(x,curU,N,dt,eps,h,S);
            [x,flag,res,iter]=pcg(func,b,tol1,400);
            %[x,flag]=minres(func,b,1e-6,100,[],[],curU(:));
            curU=curU+reshape(x, [N,N]);
            if (max(abs(curU(:)-lastU(:))))<tol
                break
            end
        end
        Uall(step+1,:,:)=curU;
    end
    U=curU;
    
end
toc;
save(['results/multistepac2s' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
end

else
    
[gamma,m]=getgamma(3,0);
maxiter=200;

tol=1e-9;
tol1=1e-3;
T=20;
eps=1;

nts=[2^3,2^4,2^5,2^6,2^7,2^8,2^9,2^10];


for numt=1:9
    %N=ceil(512*2^((numt-2)/2));
    N=2048;
    Uall=zeros(m+1,N,N);
    [Uinit,X,Y,h]=initializebigcircle(N,eps);
    %force=1/4*1./(1+exp(-1/eps^2*(4-sqrt(X.^2+Y.^2))));
    force=0;
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
            b=constructb(curU,Un,N,dt,eps,h,S,force);
            func=@(x) applyoperator(x,curU,N,dt,eps,h,S);
            [x,flag,res,iter]=pcg(func,b,tol1,400);
            %[x,flag]=minres(func,b,1e-6,100,[],[],curU(:));
            curU=curU+reshape(x, [N,N]);
            if (max(abs(curU(:)-lastU(:))))<tol
                break
            end
        end
        Uall(step+1,:,:)=curU;
    end
    U=curU;
    
end
toc;
save(['results/multistepac3s' num2str(nt) 'N' num2str(N)],'T','dt','N','U','eps');
end

end