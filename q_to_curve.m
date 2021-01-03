function p = q_to_curve(q)

[n,T] = size(q);

for i = 1:T
    qnorm(i) = norm(q(:,i),'fro');
end

for i = 1:n
    p(i,:) = [ cumtrapz(linspace(0,1,T), q(i,:).*qnorm ) ] ;
end