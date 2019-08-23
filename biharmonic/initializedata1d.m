function [uinit,ut,x,h]=initializedata1d(N,T)
h=2*pi/(N);
x=-pi:h:pi-h;

uinit=cos(x);

ut=real(1/2*exp(-T-1j*x).*(1+exp(2*1j*x)));