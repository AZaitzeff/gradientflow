addpath('../../../MATLAB/PlotPub/lib/')

N=2048;
T=1;
[Uinit,trueU,x,h]=initializedata1d(N,T);


plot(x,trueU);hold on;
plot(x,Uinit);
plt=Plot();
plt.Colors = {[.5,.5,.5],[0,0,0]};
plt.XLabel='x';
plt.YLabel='y';
plt.XLim=[-pi,pi];

plt.export('Biharmonic.png');