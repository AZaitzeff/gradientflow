function val=endfunc(x,gamma)

gamma(end,:)=x;

[cond,flag]=checkunconditionallystability(gamma);
if flag
    val=sum(-log(cond));
else
    val=inf;
end