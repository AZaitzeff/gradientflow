function gamma=xtogamma(x,m)
% gamma=zeros(m,m);
% for i=1:m
%     start=i*(i-1)/2+1;
%     gamma(i,1:i)=x(start:start+i-1);
% end

gamma=zeros(m,m);
gamma(1,1)=1;
for i=2:m
 start=(i-1)*(i-2)/2+1;
 gamma(i,2:i)=x(start:start+i-2);
 gamma(i,1)=1-sum(x(start:start+i-2));
end