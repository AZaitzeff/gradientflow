

%%
%E(u)=cosh(u)
addpath('odetests/')
order=2;

%error1 is Table 1 (order=2) or Table 2 (order=3)  
cothtest

%%
%non-smooth

addpath('odetests/')
order=2;

%error1 is the figure 3 (order=2) or Table 4 (order=3)  
quadratictest

%%
%Heat equation
addpath('heatequation/')
order=2;

%error1 is Table 3 (order=2) or Table 4 (order=3)  
mlstepfft1d



%%
%1d Allen-Cahn
addpath('allencahn/')

order=3;


%error1 is Table 5 (order=2) or Table 6 (order=3) 
multisteponedac


%%
%2d Allen-Cahn (takes a long time to run)
addpath('allencahn/')

order=2;
sieac %Creates 'true' solution
mlac

%error1 is Table 7 (order=2) or Table 8 (order=3) 
createerrortablemlac





%%
%Cahn-Hillard (takes a really long time to run)
addpath('cahnhilliard/')

order=2;
siech %Creates 'true' solution
mlch

%error1 is Table 9 (order=2) or Table 10 (order=3) 
createerrortablemlch
