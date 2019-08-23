function [A,B,C,D,E,F,G]=getconsistvars(gamma,theend)
N=size(gamma,1);
A=zeros(1,N);
B=zeros(1,N);
C=zeros(1,N);
D=zeros(1,N);
E=zeros(1,N);
F=zeros(1,N);
G=zeros(1,N);
A(1)=1/gamma(1,1);
B(1)=1/gamma(1,1)^2;
C(1)=1/2*1/gamma(1,1)^3;
D(1)=1/gamma(1,1)^3;
E(1)=1/6*1/gamma(1,1)^4;
F(1)=3/2*1/gamma(1,1)^4;
G(1)=1/gamma(1,1)^4;
for i=2:N
    S=sum(gamma(i,1:N));
    theA=sum(A(1:i-1).*gamma(i,2:i));
    theB=sum(B(1:i-1).*gamma(i,2:i));
    theC=sum(C(1:i-1).*gamma(i,2:i));
    theD=sum(D(1:i-1).*gamma(i,2:i));
    theE=sum(E(1:i-1).*gamma(i,2:i));
    theF=sum(F(1:i-1).*gamma(i,2:i));
    theG=sum(G(1:i-1).*gamma(i,2:i));
    A(i)=(1+theA)/S;
    B(i)=(A(i)+theB)/S;
    C(i)=(A(i)^2/2+theC)/S;
    D(i)=(B(i)+theD)/S;
    E(i)=(A(i)^3/6+theE)/S;
    F(i)=(A(i)*B(i)+C(i)+theF)/S;
    G(i)=(D(i)+theG)/S;
    
end

if theend==1
    A=A(N);
    B=B(N);
    C=C(N);
    D=D(N);
    E=E(N);
    F=F(N);
    G=G(N);
end

end