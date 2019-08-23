function [uinit,ut,X,Y,h]=initializedata2d(N,T)
h=40/(N);
x=-20:h:20-h;
[X,Y]=meshgrid(x);
uinit=exp(-(X.^2+Y.^2)/4);

ut=1/sqrt(T+1)*exp(-(X.^2+Y.^2)/(4+4*T));