function [gam0I,energy,iter] = InitGam(I,n,N)
% Inputs:
% I = landmark location parameter value indices
% n = half the # of basis elements to construct
% N = length of q vectors
%
% Outputs:
% gam0 = initial re-parameterization to use for main algorithm
% energy, iter

% Find landmark location parameter values
t = linspace(0,1,N);
s(1,:) = t(I{1});
s(2,:) = t(I{2});

% 
% gamint = linspace(0,t(I{2}(1)),I{1}(1));
% gamint = [gamint,linspace(t(I{2}(1)+1),t(I{2}(2)),I{1}(2)-I{1}(1))];
% gamint = [gamint,linspace(t(I{2}(2)+1),t(I{2}(3)),I{1}(3)-I{1}(2))];
% gamint = [gamint,linspace(t(I{2}(3)+1),1,N-I{1}(3))];

% Build gamint - linear interpolation through the landmark constraints
gamint = linspace(0,t(I{2}(1)),I{1}(1));

for i=1:(size(I{2},1)-1)
    gamint = [gamint,linspace(t(I{2}(i)+1),t(I{2}(i+1)),I{1}(i+1)-I{1}(i))];
end

gamint = [gamint,linspace(t(I{2}(end)+1),1,N-I{1}(end))];

% Initialize search for initial gamma with identity parameterization
iter = 1;
gam(iter,:) = linspace(0,1,N);

% Generate basis
B = IdentityBasis(n,N);
dimB = size(B,2);

% Gradient-descent algorithm
eps = 0.1;    % step size
tol = 1e-5;   % cut-off tolerance for energy

lambda = 1;

% Initialize
% No penalty
%E = norm(gam(iter,I{2})-s(1,:))^2;
%energy(iter) = E;

% Derivative penalty
% E = (norm(gam(iter,I{1})-s(2,:))^2)+lambda*trapz(linspace(0,1,N),gradient(gam(iter,:),1/(N-1)).*gradient(gam(iter,:),1/(N-1)));
% energy(iter) = E;

% Distance to gamint penalty
E = (norm(gam(iter,I{1})-s(2,:))^2)+lambda*(trapz(linspace(0,1,N),gam(iter,:)-gamint))^2;
energy(iter) = E;

if E < tol
    gam0I = gam(1,:);
else
    while E > tol
        gradE = 0;
        iter = iter+1;
        
        % Inner product for gradient descent
        % No penalty
        %innprod = dot(B(I{1},:),repmat((gam(iter-1,I{1})-s(2,:))',1,dimB));
        
        % Derivative penalty
        %innprod = dot(B(I{1},:),repmat((gam(iter-1,I{1})-s(2,:))',1,dimB))+lambda*trapz(linspace(0,1,N),(repmat(gradient(gam(iter-1,:),1/(N-1)),dimB,1).*gradient(B',1/(N-1)))');
        
        % Distance to gamint penalty
        innprod = dot(B(I{1},:),repmat((gam(iter-1,I{1})-s(2,:))',1,dimB))+lambda*trapz(linspace(0,1,N),(repmat(gam(iter-1,:)-gamint,dimB,1).*B')');
        
        gradE = gradE + sum(repmat(innprod,N,1).*B,2);
        gam(iter,:) = gam(iter-1,:)-eps*gradE';

        % No penalty
        %E = norm(gam(iter,I{2})-s(1,:))^2;
        %energy(iter) = E;
        
        % Derivative penalty
        % E = (norm(gam(iter,I{1})-s(2,:))^2)+lambda*trapz(linspace(0,1,N),gradient(gam(iter,:),1/(N-1)).*gradient(gam(iter,:),1/(N-1)));
        % energy(iter) = E;
        
        % Distance to gamint penalty
        E = (norm(gam(iter,I{1})-s(2,:))^2)+lambda*(trapz(linspace(0,1,N),gam(iter,:)-gamint))^2;
        energy(iter) = E;
        
%         plot(gam(iter,:))
%         energy(iter)
%         pause
    end
    gam0I = gam(end,:);
end
