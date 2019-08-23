function x=gammatox2(gamma,m)
x=zeros(2*m-1);
x(1)=gamma(1,1);
for i=2:m
    x(2*i-2)=gamma(i,1);
    x(2*i-1)=gamma(i,i);
end