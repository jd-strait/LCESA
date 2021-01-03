function val = InnerProd_Basis(x1,x2)

[n,N,t] = size(x1);

% val = trapz(linspace(0,1,t),trapz(linspace(0,1,N),sum(x1.*x2),2)));
val = sum(trapz(linspace(0,1,N),sum(x1.*x2),2))/(t-1);