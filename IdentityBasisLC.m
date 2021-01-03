function [B,Bnew,BnewO] = IdentityBasisLC(n,N,I)
% Inputs:
% n = # of basis elements
% N = length of q vector
% I = landmark location parameter values of curve
%
% Outputs:
% B = matrix of original basis elements
% Bnew = matrix of landmark-constrained basis elements (not orthonormal)
% BnewO = matrix of orthonormal landmark-constrained basis elements

% Get initial elements used to construct new basis
k = length(I);
B = IdentityBasis(floor(n/2),N);

% Find landmark location parameter values
t = linspace(0,1,N);
s = t(I);

% Create new basis which is zero at landmark locations
Bnew = B;

% for i=1:k
%     n = size(Bnew,2);
%     keyboard;
%     for j=1:(n-1)
%         Bnew(:,j) = Bnew(:,j)-(Bnew(I(i),j)/Bnew(I(i),end))*Bnew(:,end);
%     end
%     Bnew(:,end) = [];
% end

tol = 1e-3;
idx = 1;

for i=1:k
    n = size(Bnew,2);
    
    for m=1:n
        if abs(Bnew(I(i),m)) > tol
            idx = min(idx,m);
        else
            idx = idx + 1;
        end
    end
    
    BnewE = Bnew(:,idx);
%     for j=2:n
%         Bnew(:,j) = Bnew(:,j)-(Bnew(I(i),j)/Bnew(I(i),1))*Bnew(:,1);
%     end

    for j=1:n
        Bnew(:,j) = Bnew(:,j)-(Bnew(I(i),j)/BnewE(I(i)))*BnewE;
    end
    
    Bnew(:,idx) = [];
end

% Make the basis orthonormal via Gram-Schmidt under Palais metric
n = size(Bnew,2);
BnewO(:,1) = Bnew(:,1);

for i=2:n
    proj = 0;
    for j=1:(i-1)
    	proj = proj + (PalaisMetric(Bnew(:,i),BnewO(:,j))/PalaisMetric(BnewO(:,j),BnewO(:,j)))*BnewO(:,j);
    end
    BnewO(:,i) = Bnew(:,i)-proj;
end

for i=1:n
    BnewO(:,i) = BnewO(:,i)/sqrt(PalaisMetric(BnewO(:,i),BnewO(:,i)));
end

% keyboard;

% Check that new orthonormal basis goes through zero at landmark locations
% for i=1:n
%     clf
%     hold on
%     plot(t,BnewO(:,i))
%     scatter(t(I),BnewO(I,i))
%     pause
% end

BnewO(I,:) = 0;