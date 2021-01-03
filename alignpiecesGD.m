function [X1n,q1n,X2n,q2n,gam,energy1,energy2] = alignpiecesGD(X1,X2,I,N)

% Find q functions and optimal rotation
q1 = curve_to_q(X1);
q2 = curve_to_q(X2);

q1 = q1/sqrt(InnerProd_Q(q1,q1));
q2 = q2/sqrt(InnerProd_Q(q2,q2));

A = q1*q2';
[U,~,V] = svd(A);
if det(A)> 0
    Ot = U*V';
else
    if (size(X1,1)==2)
        Ot = U*([V(:,1) -V(:,2)])';
    else
        Ot = U*([V(:,1) V(:,2) -V(:,3)])';
    end
end

q2n = Ot*q2;
X2n = Ot*X2;

% Find optimal re-parameterization
[~,gam,energy1,energy2] = GDReparLC(X1,X2n,N,I);

q2n = Group_Action_by_Gamma_Coord_q(q2n,gam);
X2n = Group_Action_by_Gamma_Coord(X2n,gam);

% Find optimal rotation
A = q1*q2n';
[U,~,V] = svd(A);
if det(A)> 0
    Ot = U*V';
else
    if (size(X1,1)==2)
        Ot = U*([V(:,1) -V(:,2)])';
    else
        Ot = U*([V(:,1) V(:,2) -V(:,3)])';
    end
end

q2n = Ot*q2n;
X2n = Ot*X2n;

X1n = X1;
q1n = q1;