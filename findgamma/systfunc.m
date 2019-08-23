function val=systfunc(x,m,a,order,quad)

gamma=xtogamma(x,m);


[A,B,C,D,E,F,G]=getconsistvars(gamma,1);

[cond,flag]=checkunconditionallystability(gamma);
if flag
    if order==2
        val=a*((B/A^2-1/2)^2)+sum(-log(cond));
    elseif order==3
        if quad==1
            val=a*((A-1)^2+(B-1/2)^2+(D-1/6)^2)+sum(-log(cond));
        else
            val=a*((A-1)^2+(B-1/2)^2+(C-1/6)^2+(D-1/6)^2)+sum(-log(cond));
        end
    elseif order==4
        if quad==1
            val=a*((A-1)^2+(B-1/2)^2+(D-1/6)^2+(G-1/24)^2)+sum(-log(cond));
        else
            val=a*((A-1)^2+(B-1/2)^2+(C-1/6)^2+(D-1/6)^2+(E-1/24)^2+(F-1/6)^2+(G-1/24)^2)+sum(-log(cond));
        end
    end
else
    val=inf;
end