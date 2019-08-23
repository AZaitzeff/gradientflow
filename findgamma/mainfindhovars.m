maxiter1=1;
maxiter2=10000;
m=3;
order=2;
quad=0;
warning('off','MATLAB:nearlySingularMatrix')
for zall=1:maxiter1
for z=1:maxiter2
    gamma=zeros(m,m);
    gamma(1,1)=3+randn();
    gamma(1,1)=1;
    for i=2:m
        
        gamma(i,2:i)=(randn(1,i-1));
        gamma(i,1)=1-sum(gamma(i,2:i));
        %gamma(i,1:i)=sqrt(i)*randn(1,i);
        %gamma(i,i)=5+i/2+abs(sqrt(i)*randn());
    end
    [cond,flag]=checkunconditionallystability(gamma);
    if flag
        break
    end
end

options = optimset('MaxIter',1000000,'MaxFunEvals',1000000,'Display','off');


x0=gammatox(gamma,m);
func=@(x) systfunc(x,m,1e4,order,quad);
[x,fval,exitflag] = fminsearch(func,x0,options);
for a=[1e5,1e6,1e7,1e8]
    func=@(x) systfunc(x,m,a,order,quad);
    [x,fval,exitflag] = fminsearch(func,x,options);
end
gamma=xtogamma(x,m);
% 
% [A,B,C,D,E,F,G]=getconsistvars(gamma(1:m-1,1:m-1),0);
% 
% 
% coef=zeros(7,m);
% coef(1,1)=1;
% coef(1,2:m)=1-A;
% coef(2,1)=1/2;
% coef(2,2:m)=1/2-B;
% coef(3,1)=1/6;
% coef(3,2:m)=1/6-C;
% coef(4,1)=1/6;
% coef(4,2:m)=1/6-D;
% coef(5,1)=1/24;
% coef(5,2:m)=1/24-E;
% coef(6,1)=1/6;
% coef(6,2:m)=1/6-F;
% coef(7,1)=1/24;
% coef(7,2:m)=1/24-G;
% b=[1;1;1/2;1/2;1/6;1/2+1/6;1/6];
% 
%     if order==2
%         coef=coef([1,2],:);
%         b=b([1,2],:);
%     elseif order==3
%         if quad==1
%             coef=coef([1,2,4],:);
%             b=b([1,2,4],:);
%         else
%             coef=coef([1,2,3,4],:);
%             b=b([1,2,3,4],:);
%         end
%     elseif order==4
%         if quad==1
%             coef=coef([1,2,4,7],:);
%             b=b([1,2,4,7],:);
%         end
%     end
% 
% options = optimoptions('lsqlin','Display','off');
% saveg=gamma(m,:);
% x2=lsqminnorm(coef,b);
% %x2 = lsqlin(eye(m),saveg,[],[],coef,b,[],[],[],options);
% gamma(m,:)=x2;
% 
% 
% [cond,flag]=checkunconditionallystability(gamma);
% if flag && all(abs(gamma(:))<1e5)
%     'least squares'
%     break
% end
% if all(abs(gamma(:))<1e5)
% func=@(x) endfunc(x,gamma);
% options = optimoptions('fmincon','Display','off');
% [x,fval,exitflag] =fmincon(func,saveg,[],[],coef,b,[],[],[],options);
% gamma(m,:)=x;
% 
% if exitflag>0 && fval<inf
%     'fmincon'
%     break
% end
% end
[cond,flag]=checkunconditionallystability(gamma);
[A,B,C,D,E,F,G]=getconsistvars(gamma,1);
    if order==2
        val=((A-1)^2+(B-1/2)^2);
    elseif order==3
        if quad==1
            val=((A-1)^2+(B-1/2)^2+(D-1/6)^2);
        else
            val=((A-1)^2+(B-1/2)^2+(C-1/6)^2+(D-1/6)^2);
        end
    elseif order==4
        if quad==1
            val=((A-1)^2+(B-1/2)^2+(D-1/6)^2+(G-1/24)^2);
        else
            val=((A-1)^2+(B-1/2)^2+(C-1/6)^2+(D-1/6)^2+(E-1/24)^2+(F-1/6)^2+(G-1/24)^2);
        end
    end

if val<1e-2 && flag && all(abs(gamma(:))<1000)
   break
end
end
%save(['found' num2str(m)],'gamma');
warning('on','MATLAB:nearlySingularMatrix')