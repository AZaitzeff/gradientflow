n=1001;
h=4/(n-1);
x=-2:h:2;
eps=1/4;
c=2;
T=1/2;
uinit=tanh((x+.5)/eps);
trueu=tanh(((x+.5)-c*T)/eps);
nts=[2^3,2^4,2^5,2^6,2^7,2^8]*6;
error1=zeros(2,6);

tol=1e-10;
%U=zeros(m+1,n);
maxiter=100;
for numt=1:6
    nt=nts(numt);
    dt=T/nt;
    u=uinit;
    for t=1:nt
    utemp=u;
    ucur=u;
    for z =1:maxiter
        ulast=ucur;
        B=zeros(n-4,5);
        B(:,3)=5/2*dt/h^2+1+dt*(-2/eps^2-2*c/eps*ucur(3:n-2)+6/eps^2*ucur(3:n-2).^2);
        B(:,2)=-4/3*dt/h^2;
        B(:,4)=-4/3*dt/h^2;
        B(:,1)=1/12*dt/h^2;
        B(:,5)=1/12*dt/h^2;
        F=-dt/h^2*(-5/2*ucur(3:n-2)+4/3*(ucur(2:n-3)+ucur(4:n-1))-1/12*(ucur(1:n-4)+ucur(5:n)))+ucur(3:n-2)...
            +dt*(c/eps-2*ucur(3:n-2)/eps^2-c/eps*ucur(3:n-2).^2+2/eps^2*ucur(3:n-2).^3)-utemp(3:n-2);
        A = spdiags(B,[-2,-1,0,1,2],n-4,n-4);
        ucur(3:n-2) = ucur(3:n-2)- mldivide(A,F')';

        if sum(abs(ulast-ucur))<tol
            break
        end
    end
    u=ucur;
    end
    error1(1,numt)=sqrt(sum((trueu-u).^2))/sqrt(sum((trueu-uinit).^2));
end

for i=2:6
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end