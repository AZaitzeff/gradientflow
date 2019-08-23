function [cond,flag]=checkunconditionallystability(gamma)
S=sum(gamma,2);
n=numel(S);
cond=zeros(n,1);


for i=n:-1:1
   S(i)=S(i)-sum(S(i+1:n));
   cond(i)=S(i);
   for j=1:i
      gamma(i,j)=(gamma(i,j)-sum(S(i+1:n).*gamma(i+1:n,j)))/S(i); 
   end
   for z=i:n
        gamma(z,i)=1-gamma(z,i);
        S(z)=S(z)*(gamma(z,i))^2;
        for j=1:i-1
            gamma(z,j)=gamma(z,j)/(gamma(z,i));
        end
        gamma(z,i)=0;
   end
end

if any(cond<0)
    flag=0;
else
    flag=1;
end