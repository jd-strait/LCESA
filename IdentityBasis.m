function B = IdentityBasis(n,N)
% Inputs:
% n = half the # of basis elements desired
% N = length of q vector
% Output:
% B = matrix of basis elements

t = linspace(0,1,N);

for i=1:n
    B(:,i) = (1/(sqrt(2)*pi*i))*sin(2*pi*i*t);
    B(:,n+i) = (1/(sqrt(2)*pi*i))*(cos(2*pi*i*t)-1);
end