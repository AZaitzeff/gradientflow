function [U,h]=initializecircle(N,eps)
h=2/(N);
x=-1:h:1-h;
[X,Y] = meshgrid(x);
U=1./(1+exp(-1/eps*(.5-sqrt(X.^2+Y.^2))));