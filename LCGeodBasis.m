function bO = LCGeodBasis(n1,N1,n2,N2)
%% Inputs
% n1 = half the number of periodic basis elements for x, y
% N1 = number of points sampled on each shape
% n2 = half the number of basis elements for t
% N2 = number of points along geodesic path to sample

%% Output
% b0 = orthonormal basis for landmark-constrained geodesic (for path
% straightening algorithm)

%% Spatial basis
n1 = 10; % half the number of periodic basis elements for x, y
N1 = 100; % number of points to sample for curve

t1 = linspace(0,1,N1);

% Form basis in x dimension
for i=1:n1
    Bx(1,:,i) = t1;
    Bx(2,:,i) = sin(2*pi*i*t1);
    Bx(1,:,n1+i) = t1;
    Bx(2,:,n1+i) = cos(2*pi*i*t1);
    Bx(1,:,2*n1+1) = t1;
    Bx(2,:,2*n1+1) = t1;
    Bx(1,:,2*n1+2) = t1;
    Bx(2,:,2*n1+2) = 1-t1;
end

for i=1:(2*n1+2)
    Bx(3,:,i) = zeros(1,N1);
end

% Form basis in y dimension
By(1,:,:) = Bx(1,:,:);
By(2,:,:) = Bx(3,:,:);
By(3,:,:) = Bx(2,:,:);

% Temporal basis
n2 = 2; % half the number of basis elements for t
N2 = 7; % number of points along geodesic path to sample

t2 = linspace(0,1,N2);

% Form basis in t dimension
for i=1:n2
    Bt(1,:,i) = t2;
    Bt(2,:,i) = sin(2*pi*i*t2);
    Bt(1,:,n2+i) = t2;
    Bt(2,:,n2+i) = cos(2*pi*i*t2)-1;
end

% Combine spatial and temporal bases
for i=1:size(Bx,3)
    for j=1:size(Bt,3)
        for k=1:N2
            b1(1,:,k,2*n2*(i-1)+j) = t1;
            b1(2:3,:,k,2*n2*(i-1)+j) = Bx(2:3,:,1) * Bt(2,k,j);
        end
    end
end

for i=1:size(By,3)
    for j=1:size(Bt,3)
        for k=1:N2
            b2(1,:,k,2*n2*(i-1)+j) = t1;
            b2(2:3,:,k,2*n2*(i-1)+j) = By(2:3,:,1) * Bt(2,k,j);
        end
    end
end

% All possible combinations of b1 and b2 (this is the full basis)
for i=1:size(b1,4)
    for j=1:size(b2,4)
        b(1,:,:,size(b1,4)*(i-1)+j) = b1(1,:,:,1);
        b(2:3,:,:,size(b1,4)*(i-1)+j) = b1(2:3,:,:,i) + b2(2:3,:,:,j);
    end
end

% Orthonormalize via Gram-Schmidt under L2 metric
%n = size(b,4);
n = 50;
bO(:,:,:,1) = b(:,:,:,1);

for i=2:n
    proj = 0;
    for j=1:(i-1)
    	proj = proj + InnerProd_3d(b(:,:,:,i),bO(:,:,:,j))./InnerProd_3d(bO(:,:,:,j),bO(:,:,:,j)).*bO(:,:,:,j);
    end
    bO(:,:,:,i) = b(:,:,:,i)-proj;
    bO(1,:,:,i) = repmat(linspace(0,1,N1),1,1,N2,1);
end

for i=1:n
    bO(:,:,:,i) = bO(:,:,:,i)./sqrt(InnerProd_3d(bO(:,:,:,i),bO(:,:,:,i)));
    bO(1,:,:,i) = repmat(linspace(0,1,N1),1,1,N2,1);
end