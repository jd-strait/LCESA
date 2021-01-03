function fn = Group_Action_by_Gamma_Coord_q(f,gamma)

[n,T] = size(f);

gam_dev = gradient(gamma, 1/T);

for j=1:n
    fn(j,:) = interp1(linspace(0,1,T) , f(j,:),gamma,'linear').*sqrt(gam_dev);
end