function[X1n,q1n,X2n,q2n,q1l,q2ln,gam] = alignDP(X1,X2,Y1,Y2)

% Load some parameters, no need to change this
    lam = 0;
    N = 100;
    
    numPieces = size(Y1,1);
    q1l = cell(numPieces,1);
    q2l = cell(numPieces,1);
    l1 = 0; l2 = 0;
    for i=1:numPieces
        Y1{i} = ReSampleCurve(Y1{i},N);
        q1l{i} = curve_to_q2(Y1{i});
        l1 = l1 + InnerProd_Q(q1l{i},q1l{i});
    end
    for j=1:numPieces
        Y2{j} = ReSampleCurve(Y2{j},N);
        q2l{j} = curve_to_q2(Y2{j});
        l2 = l2 + InnerProd_Q(q2l{j},q2l{j});
    end

    for i=1:numPieces
        q1l{i} = q1l{i}/sqrt(l1);
        q2l{i} = q2l{i}/sqrt(l2);
    end
    
% Form the q function for representing curves and find best rotation
    q1 = curve_to_q(X1);
    q2 = curve_to_q(X2);
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
    
    for j=1:size(Y2,1)
        q2l{j} = Ot*q2l{j};
        Y2{j} = Ot*Y2{j};
    end

% Realigning pieces of curves, and putting pieces back together
    gamI = zeros(numPieces,N);
    Y2n = cell(numPieces,1);
    q2ln = cell(numPieces,1);
    l2n = 0;
    for i=1:numPieces
        gamI(i,:) = DynamicProgrammingQ(q2l{i}/sqrt(InnerProd_Q(q2l{i},q2l{i})),q1l{i}/sqrt(InnerProd_Q(q1l{i},q1l{i})),lam,0);
        keyboard;
        gamI(i,:) = (gamI(i,:)-gamI(i,1))/(gamI(i,end)-gamI(i,1));
        Y2n{i} = Group_Action_by_Gamma_Coord(Y2{i},gamI(i,:));
        q2ln{i} = curve_to_q2(Y2n{i});
        l2n = l2n + InnerProd_Q(q2ln{i},q2ln{i});
    end
    

    
    for i=1:numPieces
        q2ln{i} = q2ln{i}/sqrt(l2n);
    end
   
    X2n = zeros(2,(numPieces-1)*(N-1) + N);
    X1n = zeros(2,(numPieces-1)*(N-1) + N);
    
    n = size(Y1{1},1);
    for j=1:numPieces
%         for i = 1:n
%             v1(i,:) = gradient(Y1{j}(i,:),1/(N-1));
%             v2(i,:) = gradient(Y2{j}(i,:),1/(N-1));
%         end
        len1(j) = InnerProd_Q(q1l{j},q1l{j});
        len2(j) = InnerProd_Q(q2l{j},q2l{j});
    end

    for j=1:numPieces
        gamI(j,:)=gamI(j,:)*len2(j);
    end
    
    gam(1,:)=linspace(0,len1(1),N);
    gam(2,:)=gamI(1,:);
    for j=2:numPieces
        gam(1,N*(j-1)+1:j*N)=linspace(sum(len1(1:j-1)),sum(len1(1:j)),N);
        gam(2,N*(j-1)+1:j*N)=sum(gamI(1:j-1,end))+gamI(j,:);
    end
    
    for j=1:numPieces
        if j==numPieces
            X2n(:,(N*j - N - j + 2):(((N-1)*j))+1) = Y2n{j}(:,1:N);
            X1n(:,(N*j - N - j + 2):(((N-1)*j))+1) = Y1{j}(:,1:N);
        else
            X2n(:,(N*j - N - j + 2):((N-1)*j)) = Y2n{j}(:,1:(N-1));
            X1n(:,(N*j - N - j + 2):((N-1)*j)) = Y1{j}(:,1:(N-1));
        end
    end

    q2n = curve_to_q(X2n);
    q1n = curve_to_q(X1n);

%Find optimal rotation
    A = q1n*q2n';
    [U,S,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        if (size(X1,1)==2)
            Ot = U*([V(:,1) -V(:,2)])';
        else
            Ot = U*([V(:,1) V(:,2) -V(:,3)])';
        end
    end
    X2n = Ot*X2n;
    q2n = Ot*q2n;

    for j=1:size(Y2,2)
        q2ln{j} = Ot*q2ln{j};
    end