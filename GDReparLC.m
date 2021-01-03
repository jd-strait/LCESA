function [gam0I,gamF,energy1,energy2] = GDReparLC(X1,X2,N,I)
% Inputs:
% X1, X2 = 2 curves (two-dimensional)
% N = number of points to resample each curve to
% I = indices of landmarks
%
% Outputs:
% gam0 = initialized re-parameterization
% gamF = optimal re-parameterization
% energy1, energy 2

% X1 = ReSampleCurve(X1,N);
% X2 = ReSampleCurve(X2,N);

% Convert to q (unit length)
q1 = curve_to_q(X1);
q2 = curve_to_q(X2);

% % Select landmarks - I is cell of indices of landmark locations on each
% % curve
% [~,I] = cursor({X1,X2});

% Initialize gamma
[gam0I,energy1,~] = InitGam(I,floor(N/2),N);

% Find optimal gamma in landmark-constrained re-parameterization group
[gamF,energy2,~,q2t] = OptRepLC(gam0I,I,40,q1,q2);