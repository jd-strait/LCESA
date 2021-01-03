function BstO = BasisGeod(q1,n1,n2,N2)

N1 = size(q1,2); % number of points to sample for curve

t1 = linspace(0,1,N1);

% Form basis in x dimension
Bx(1,:,1) = t1;
Bx(1,:,2) = 1-t1;

for i=1:n1
    Bx(1,:,i+2) = sin(2*pi*i*t1);
    Bx(1,:,n1+i+2) = cos(2*pi*i*t1);
end

for i=1:(2*n1+2)
    Bx(2,:,i) = zeros(1,N1);
end

% Form basis in y dimension
By(1,:,:) = Bx(2,:,:);
By(2,:,:) = Bx(1,:,:);

% Combine spatial bases
n = 1;

for i=1:size(Bx,3)
    Bs(:,:,n) = Bx(:,:,i);
    n = n+1;
end

for j=1:size(By,3)
    Bs(:,:,n) = By(:,:,j);
    n = n+1;
end

% Temporal basis
%n2 = 2; % half the number of basis elements for t
%N2 = 7; % number of points along geodesic path to sample

t2 = linspace(0,1,N2);

% Form basis in t dimension
for i=1:n2
    Bt(1,:,i) = sin(2*pi*i*t2);
    Bt(1,:,n2+i) = cos(2*pi*i*t2)-1;
end

% Remove basis elements that are approximately zero
idx = [];
for i=1:size(Bt,3)
    if sum(abs(Bt(:,:,i))<1e-10)==size(Bt,2)
        idx = [idx,i];
    end
end

Bt(:,:,idx) = [];

% Zero out the final component to avoid numerical issues
for i=1:size(Bt,3)
    Bt(1,N2,i) = 0;
end

% Combine spatial basis with temporal basis
for i=1:size(Bs,3)
    for j=1:size(Bt,3)
        for k=1:N2
            Bst(:,:,k,size(Bt,3)*(i-1)+j) = Bs(:,:,i) * Bt(:,k,j);
        end
    end
end

% Orthonormalize via Gram-Schmidt under L2 metric
n = size(Bst,4);
%n = 400;
BstO(:,:,:,1) = Bst(:,:,:,1);

for i=2:n
    proj = zeros(2,N1,N2);
    for j=1:(i-1)
        if abs(InnerProd_Basis(BstO(:,:,:,j),BstO(:,:,:,j))) > 1e-16
            proj = proj + (InnerProd_Basis(Bst(:,:,:,i),BstO(:,:,:,j))/InnerProd_Basis(BstO(:,:,:,j),BstO(:,:,:,j)))*BstO(:,:,:,j);
        end
    end
    BstO(:,:,:,i) = Bst(:,:,:,i)-proj;
end

for i=1:n
    BstO(:,:,:,i) = BstO(:,:,:,i)/sqrt(InnerProd_Basis(BstO(:,:,:,i),BstO(:,:,:,i)));
end