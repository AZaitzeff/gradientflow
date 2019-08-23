function x=gammatox(gamma,m)
m1=m-1;
x=zeros(1,m1*(m1+1)/2);
for i=2:m
    start=(i-1)*(i-2)/2+1;
    x(start:start+i-2)=gamma(i,2:i);
end