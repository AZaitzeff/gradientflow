function [A,B,C,D]=getconsistvars(gamma,theend)
N=size(gamma,1);
A=zeros(1,N);
B=zeros(1,N);
C=zeros(1,N);
D=zeros(1,N);

A(1)=1/gamma(1,1);
B(1)=1/gamma(1,1)^2;
C(1)=1/2*1/gamma(1,1)^3;
D(1)=1/gamma(1,1)^3;

for i=2:N
    S=sum(gamma(i,1:N));
    theA=sum(A(1:i-1).*gamma(i,2:i));
    theB=sum(B(1:i-1).*gamma(i,2:i));
    theC=sum(C(1:i-1).*gamma(i,2:i));
    theD=sum(D(1:i-1).*gamma(i,2:i));

    A(i)=(1+theA)/S;
    B(i)=(A(i)+theB)/S;
    C(i)=(A(i)^2/2+theC)/S;
    D(i)=(B(i)+theD)/S;

    
end

if theend==1
    A=A(N);
    B=B(N);
    C=C(N);
    D=D(N);
end

end