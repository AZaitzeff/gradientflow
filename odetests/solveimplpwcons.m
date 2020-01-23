function un1=solveimplpwcons(k,un)

val1=un-k;
val2=un-k*1/2;


if val1>=1
    un1=val1;
elseif val1<1 && val2>1
    un1=1;
elseif val2<=1 && val2>=0
    un1=val2;
elseif val2<0
    un1=0;
else
    'error'
end