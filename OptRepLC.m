function [gamF,energy,iter,q2t] = OptRepLC(gam0,I,n,q1,q2)
% Inputs:
% gam0 = initial guess from initialization step
% I = landmark location parameter value indices
% n = # of basis elements to construct
% q1,q2 = SRVFs of 2 curves
%
% Outputs:
% gamF = optimal re-parameterization
% energy, iter

N = length(q1);

% Generate basis
[~,~,BO] = IdentityBasisLC(n,N,I{2});

% dB = size(BO,2);
% BO1 = BO(:,1:floor(dB/2));
% BO2 = BO(:,1:dB);

% Initial search
% Gradient-descent algorithm
eps = 0.005;    % step size
tol = 0.01;    % cut-off tolerance for gradient
maxIter = 15000;

% Initialize search for gamma with initialized re-parameterization
iter = 1;

gam(iter,:) = gam0;

q2t = Group_Action_by_Gamma_Coord_q(q2,gam(iter,:));
energy(iter) = InnerProd_Q(q1-q2t,q1-q2t);

dimB = size(BO,2);
d = size(q1,1);

gradH = zeros(1,N);
for i=1:dimB
    gradH = gradH + InnerProd_Q(q1-q2t,gradient(q2t,1/N).*repmat(BO(:,i)',d,1)+0.5*q2t.*repmat(gradient(BO(:,i),1/N)',d,1))*BO(:,i)';
end

ngradH = norm(gradH);

if ngradH < tol
    gamF = gam(1,:);
else
    while ngradH > tol && iter <= maxIter
        iter = iter + 1;
        inc_gam = linspace(0,1,N)+eps*gradH;
        gam(iter,:) = interp1(linspace(0,1,N),gam(iter-1,:),inc_gam);

        %plot(gam(iter,:))
        %pause
        q2t = Group_Action_by_Gamma_Coord_q(q2,gam(iter,:));
        energy(iter) = InnerProd_Q(q1-q2t,q1-q2t);
        gradH = zeros(1,N);
        for i=1:dimB
            gradH = gradH + InnerProd_Q(q1-q2t,gradient(q2t,1/N).*repmat(BO(:,i)',d,1)+0.5*q2t.*repmat(gradient(BO(:,i),1/N)',d,1))*BO(:,i)';
        end
        ngradH = norm(gradH);
    end
    gamF = gam(end,:);
end

% clf
% hold on
% figure(1);
% plot(gam0)
% plot(gamF)
% 
% figure(2);
% plot(energy)

% % Refinement
% eps = 0.01;    % step size
% tol = 0.0001;    % cut-off tolerance for gradient
% maxIter = 8000;
% 
% % Initialize search for gamma with initialized re-parameterization
% iter2 = 1;
% 
% gam2(iter2,:) = gamF;
% 
% q2t = Group_Action_by_Gamma_Coord_q(q2,gam2(iter2,:));
% energy2(iter2) = InnerProd_Q(q1-q2t,q1-q2t);
% 
% dimB = size(BO2,2);
% d = size(q1,1);
% 
% gradH = zeros(1,N);
% for i=1:dimB
%     gradH = gradH + InnerProd_Q(q1-q2t,gradient(q2t,1/N).*repmat(BO2(:,i)',d,1)+0.5*q2t.*repmat(gradient(BO2(:,i),1/N)',d,1))*BO2(:,i)';
% end
% 
% ngradH = norm(gradH);
% 
% if ngradH < tol
%     gamF1 = gam2(1,:);
% else
%     while ngradH > tol && iter2 <= maxIter
%         iter2 = iter2 + 1;
%         inc_gam = linspace(0,1,N)+eps*gradH;
%         gam2(iter2,:) = interp1(linspace(0,1,N),gam2(iter2-1,:),inc_gam);
% 
%         %plot(gam(iter,:))
%         %pause
%         q2t = Group_Action_by_Gamma_Coord_q(q2,gam2(iter2,:));
%         energy2(iter2) = InnerProd_Q(q1-q2t,q1-q2t);
%         gradH = zeros(1,N);
%         for i=1:dimB
%             gradH = gradH + InnerProd_Q(q1-q2t,gradient(q2t,1/N).*repmat(BO2(:,i)',d,1)+0.5*q2t.*repmat(gradient(BO2(:,i),1/N)',d,1))*BO2(:,i)';
%         end
%         ngradH = norm(gradH);
%     end
%     gamF1 = gam2(end,:);
% end
% 
% keyboard;
% 
% energy = [energy,energy2];