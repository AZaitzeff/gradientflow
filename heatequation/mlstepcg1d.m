tol=1e-8;
tol1=1e-4;
construct3rdordergamma;

N=2048;
T=1;
[Uinit,trueU,x,h]=initializedata1d(N,T);


l=(0:(N-1));

nts=[2,4,8,16,32,64,128,256];
NT=numel(nts);
error1=zeros(2,NT);
error2=zeros(2,NT);
error3=zeros(2,NT);


m=6;
Uall=zeros(m+1,N);
for numt=1:NT
    U=Uinit;
    nt=nts(numt);
    dt=T/nt;
    
    tic;
for t=1:nt
    curU=U;
    Uall(1,:,:)=U;
    
    for step=1:m
        Un=squeeze(Uall(1,:)*gamma(step,1));
        for zm=2:step
            Un=Un+squeeze(Uall(zm,:,:)*gamma(step,zm));
        end
        S=sum(gamma(step,:));
        
        for i=1:maxiter
            lastU=curU;
            b=-(-Un'+S*curU'-dt*laplacian51d(curU,N,h)');
            func=@(x) (S*x-dt*laplacian51d(x',N,h)');
            [x,flag,res,iter]=pcg(func,b,tol1,200);
            %[x,flag]=minres(func,b,1e-6,100,[],[],curU(:));
            curU=curU+x';
            if (max(abs(curU(:)-lastU(:))))<tol
                break
            end
        end
        
        %val=normdiff(curU,Un,N,dt,eps,h,S);
        Uall(step+1,:,:)=curU;
    end
    U=curU;
    
    
    
end
toc;
    error1(1,numt)=sqrt(h*sum((U(:)-trueU(:)).^2));
    error2(1,numt)=sqrt(h*sum((U(:)-trueU(:)).^2))/sqrt(h*sum((trueU(:)).^2));
    error3(1,numt)=max(abs(U(:)-trueU(:)));

end

for i=2:NT
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error2(2,i)=(log(error2(1,i-1))-log(error2(1,i)))/log(2);
    error3(2,i)=(log(error3(1,i-1))-log(error3(1,i)))/log(2);
end