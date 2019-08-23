function gamma=xtogamma2(x,m)
gamma=zeros(m,m);
gamma(1,1)=x(1);
for i=2:m
    gamma(i,1)=x(2*i-2);
    gamma(i,i)=x(2*i-1);
end