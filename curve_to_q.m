function [q] = curve_to_q(p)

[n,N] = size(p);
for i = 1:n
    v(i,:) = gradient(p(i,:),1/(N-1));
end

% len = trapz(linspace(0,1,N),sum(v.^2));
% v = v/len;

% Calculates q
for i = 1:N
    L(i) = sqrt(norm(v(:,i),'fro'));
    if L(i) > 0.0001
        q(:,i) = v(:,i)/L(i);
    else
        q(:,i) = v(:,i)*0.0001;
    end
end

% Standardizes q to have length 1
q = q/sqrt(InnerProd_Q(q,q));