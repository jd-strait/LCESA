function p = PalaisMetric(v1,v2)
% Inputs:
% v1, v2 = two vectors
% Output:
% p = Palais Metric

N = length(v1);

prod = gradient(v1,1/N).*gradient(v2,1/N);
p = trapz(linspace(0,1,N),prod);