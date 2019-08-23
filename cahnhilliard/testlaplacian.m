 Ns=[128,256,512,1024,2048,4096,4096*2];
 error1=zeros(2,7);
 for i=1:6
     N=Ns(i);
     h=20/N;
     x=-10:h:(10-h);

[X,Y]=meshgrid(x);
u=cos(2*pi*X).*cos(2*pi*Y);
%D=builddlap(N,h);
%L=D*u(:);

L=laplacian9(u,N,h);
L=laplacian9(L,N,h);
true=4*(2*pi)^4*cos(2*pi*X).*cos(2*pi*Y);
%newtrue=true+h^2/12*4*(2*pi)^4*cos(2*pi*X).*cos(2*pi*Y);
error1(1,i)=sqrt(sum((L(:)-true(:)).^2)*h^2);
if i>1
    error1(2,i)=(log(error1(1,i-1))-log(error1(1,i)))/log(2);
end
 end