addpath('findgamma/');
tol=1e-12;
%Checks if various gammas satifies theorem 2.1 and the consistency equation
[gamma,m]=getgamma(2,0);
[S]=checkunconditionallystability2(gamma);
[A,B,~,~]=getconsistvars(gamma,1);

all(diag(S)>0) && abs(A-1)<tol && abs(B-1/2)<tol


[gamma,m]=getgamma(2,1);

[S]=checkunconditionallystability2(gamma);

[A,B,~,~]=getconsistvars(gamma,1);

all(diag(S)>0) && abs(A-1)<tol && abs(B-1/2)<tol

[gamma,m]=getgamma(3,0);

[S]=checkunconditionallystability2(gamma);

[A,B,C,D]=getconsistvars(gamma,1);

all(diag(S)>0) && abs(A-1)<tol && abs(B-1/2)<tol && abs(C-1/6)<tol && abs(D-1/6)<tol