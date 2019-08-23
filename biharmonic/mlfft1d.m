addpath('findgamma/')
[gamma,m]=getgamma(order,0);

N=2048;
T=1;
[Uinit,trueU,x,h]=initializedata1d(N,T);


l=(0:(N-1));

nts=[4,8,16,32,64,128,256,512,1024,2048];
NT=numel(nts);
error1=zeros(2,NT);
error2=zeros(2,NT);
error3=zeros(2,NT);
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
        %invmatrix=S+dt/h^4*(28/3-13*cos(2*pi*l/N)+4*cos(4*pi*l/N)-1/3*cos(6*pi*l/N));
        invmatrix1=sqrt(S)-1j*sqrt(dt)/(h^2)*(-5/2+8/3*cos(2*pi*l/N)-1/6*cos(4*pi*l/N));
        invmatrix2=sqrt(S)+1j*sqrt(dt)/(h^2)*(-5/2+8/3*cos(2*pi*l/N)-1/6*cos(4*pi*l/N));
        %invmatrix=(invmatrix2.*invmatrix1);
        Fbar=fft(Un);
        %Ubar=(Fbar./(invmatrix));
        Ubar=(Fbar./(invmatrix1))./invmatrix2;
        curU=real(ifft(Ubar));
        %val=normdiff(curU,Un,N,dt,eps,h,S);
        Uall(step+1,:,:)=curU;
    end
    U=curU;
    
    
    
end
toc;
    error1(1,numt)=sqrt(h*sum((U(:)-trueU(:)).^2));
    error2(1,numt)=h*sum(abs(U(:)-trueU(:)));
    error3(1,numt)=max(abs(U(:)-trueU(:)));

end

for i=2:NT
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
    error2(2,i)=(log(error2(1,i-1))-log(error2(1,i)))/log(2);
    error3(2,i)=(log(error3(1,i-1))-log(error3(1,i)))/log(2);
end