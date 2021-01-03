function [q,len] = curve_to_q2(p)

[n,N] = size(p);
for i = 1:n
    v(i,:) = gradient(p(i,:),1/(N-1));
end

for i = 1:N
    L(i) = sqrt(norm(v(:,i),'fro'));
    if L(i) > 0.000000001
        q(:,i) = v(:,i)/L(i);
    else
        q(:,i) = v(:,i)*0.000000001;
    end
end