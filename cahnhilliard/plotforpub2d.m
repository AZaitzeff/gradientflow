N=1024;
[Uinit,X,Y,h]=initializebigoval(N,1);
load(['results/siech262144N1024.mat'])
x=X(1,:);

imagesc(x,x,Uinit); colormap 'gray';colorbar;
print('initial2dch','-dpng')
close;
imagesc(x,x,U,[0,1]); colormap 'gray';colorbar;
print('end2dch','-dpng')
close;