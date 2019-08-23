function [uinit,ut,x,h]=initializedata1d(N,T)
h=2/(N);
x=-1:h:1-h;

uinit=sin(pi*x);

ut=sin(pi*x)*exp(-pi^2*T);