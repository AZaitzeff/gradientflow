function un1=solveimplpwcons2(gam,un)

a=0;
b=2;
tol=1e-20;
f=@(x) func(x)+1/2*gam*(x-un)^2;

gr = (sqrt(5) + 1) / 2;

c = b - (b - a) / gr;
d = a + (b - a) / gr;
while abs(c - d) > tol
    if f(c) < f(d)
        b = d;
    else
        a = c;
    end

    c = b - (b - a) / gr;
    d = a + (b - a) / gr;
end
un1=(b + a)/2;


end

function val=func(x)
if abs(x)>=1
    val=1/2+abs(x-1);
else
    val=1/2*abs(x);
end
end
