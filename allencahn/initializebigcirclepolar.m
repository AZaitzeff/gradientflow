function [U,r,h]=initializebigcirclepolar(N,eps)
h=15/(N);
r=0:h:15-h;
U=1./(1+exp(-1/eps^2*(5-r)));