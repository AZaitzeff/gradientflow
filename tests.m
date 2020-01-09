%E(u)=u^2/2

addpath('odetests/')
order=2;

%error1 is Table 1 (order=2) or Table 2 (order=3)  
quadratictest

%%
%E(u)=cosh(u)
addpath('odetests/')
order=2;

%error1 is Table 3 (order=2) or Table 4 (order=3)  
cothtest


%%
%Heat equation
addpath('heatequation/')
order=2;

%error1 is Table 5 (order=2) or Table 6 (order=3)  
mlstepfft1d



%%
%Biharmonic
addpath('biharmonic/')
order=2;

%error1 is Table 7 (order=2) or Table 8 (order=3) 
mlfft1d

%%
%1d Allen-Cahn
addpath('allencahn/')

order=3;


%error1 is Table 9 (order=2) or Table 10 (order=3) 
multisteponedac


%%
%2d Allen-Cahn (takes a long time to run)
addpath('allencahn/')

order=2;
sieac %Creates 'true' solution
mlac

%error1 is Table 11 (order=2) or Table 12 (order=3) 
createerrortablemlac





%%
%Cahn-Hillard (takes a really long time to run)
addpath('cahnhilliard/')

order=2;
siech %Creates 'true' solution
mlch

%error1 is Table 13 (order=2) or Table 14 (order=3) 
createerrortablemlch
