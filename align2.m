function[X2n,q2n,dist] = align2(q1,q1l,X2,q2l)

    numPieces = size(q1l,2);
    N = 100;
    lam = 0;
    
    q2 = curve_to_q(X2);
% Form the q function for representing curves and find best rotation
    A = q1*q2';
    [U,~,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        Ot = U*([V(:,1) -V(:,2)])';
    end
    
    X2 = Ot*X2;
    q2 = Ot*q2;
    
    for j=1:size(q2l,2)
        q2l{1,j} = Ot*q2l{1,j};
    end
    
% Realigning pieces of curves, and putting pieces back together
    gam = zeros(numPieces,N);
    
    for i=1:numPieces
        gam(i,:) = DynamicProgrammingQ(q2l{1,i}/sqrt(InnerProd_Q(q2l{1,i},q2l{1,i})),q1l{1,i}/sqrt(InnerProd_Q(q1l{1,i},q1l{1,i})),lam,0);
        gam(i,:) = (gam(i,:)-gam(i,1))/(gam(i,end)-gam(i,1));    
    end
    
   
   gamma=gam(1,:);
   for j=2:numPieces-1
       gamma = [gamma,gam(j,2:end)+(j-1)];
   end
   gamma=[gamma,gam(numPieces,2:end)+(numPieces-1)]/numPieces;
   X2n = Group_Action_by_Gamma_Coord(X2,gamma);
   
%    q2ln{1,i} = Group_Action_by_Gamma_Coord_q(q2l{1,i},gam(i,:)); 
%    q2n = q2ln{1,1}(:,1:(N-1));
%    q1n = q1l{1,1}(:,1:(N-1));
%    for j=2:numPieces-1
%        q2n = [q2n , q2ln{1,j}(:,1:(N-1))];
%        q1n = [q1n , q1l{1,j}(:,1:(N-1))];
%    end
%    q2n = [q2n , q2ln{1,numPieces}(:,1:N)];
%    q1n = [q1n , q1l{1,numPieces}(:,1:N)];

% Find optimal rotation
    q2n = curve_to_q(X2n);
   

    q2n = ProjectC(q2n);
    A = q1*q2n';
    [U,~,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        Ot = U*([V(:,1) -V(:,2)])';
    end
    X2n = Ot*X2n;
    
%    q2n = q2ln{1,1}(:,1:(N-1));
%    q1n = q1l{1,1}(:,1:(N-1));
%    for j=2:numPieces-1
%        q2n = [q2n , q2ln{1,j}(:,1:(N-1))];
%        q1n = [q1n , q1l{1,j}(:,1:(N-1))];
%    end
%    q2n = [q2n , q2ln{1,numPieces}(:,1:N)];
%    q1n = [q1n , q1l{1,numPieces}(:,1:N)];
   
   q2n = Ot*q2n;

   dist = acos(InnerProd_Q(q1,q2n));