function[X1n,q1n,X2n,q2n,q1l,q2ln] = alignnopar(X1,X2,Y1,Y2)

% Load some parameters, no need to change this
    lam = 0;
    N = 100;
    
    numPieces = size(Y1,1);
    q1l = cell(numPieces,1);
    q2ln = cell(numPieces,1);
    l1 = 0; l2 = 0;
    for i=1:numPieces
        Y1{i} = ReSampleCurve(Y1{i},N);
        q1l{i} = curve_to_q2(Y1{i});
        l1 = l1 + InnerProd_Q(q1l{i},q1l{i});
    end
    for j=1:numPieces
        Y2{j} = ReSampleCurve(Y2{j},N);
        q2ln{j} = curve_to_q2(Y2{j});
        l2 = l2 + InnerProd_Q(q2ln{j},q2ln{j});
    end

    for i=1:numPieces
        q1l{i} = q1l{i}/sqrt(l1);
        q2ln{i} = q2ln{i}/sqrt(l2);
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
        q2ln{j} = Ot*q2ln{j};
        Y2{j} = Ot*Y2{j};
    end
   
    X2n = zeros(2,(numPieces-1)*(N-1) + N);
    X1n = zeros(2,(numPieces-1)*(N-1) + N);
    
    for j=1:numPieces
        if j==numPieces
            X2n(:,(N*j - N - j + 2):(((N-1)*j))+1) = Y2{j}(:,1:N);
            X1n(:,(N*j - N - j + 2):(((N-1)*j))+1) = Y1{j}(:,1:N);
        else
            X2n(:,(N*j - N - j + 2):((N-1)*j)) = Y2{j}(:,1:(N-1));
            X1n(:,(N*j - N - j + 2):((N-1)*j)) = Y1{j}(:,1:(N-1));
        end
    end

    q2n = curve_to_q(X2n);
    q1n = curve_to_q(X1n);