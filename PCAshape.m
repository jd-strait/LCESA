% PCA on the vectors in the tangent space to the average q
% Input: v = vectors from karchermean program, stored as row array; avgq =
% q representation of mean shape
function [X,u,S,len] = PCAshape(avgq,v,loading,Closed,disp)

% n = # of curves
n=size(v,2);


for i=1:n
    vs(i,:) = [v{1,i}(1,:),v{1,i}(2,:)];
end

K = cov(vs);

% Returns largest variances and their loadings
[U,S,V] = svd(K);
% [U,S,V] = svds(K,n-1);  % If want n-1 largest variances

% Convert back to two-dimensional loadings
for i=1:size(U,2)
    u(1,:,i) = U(1:(0.5*size(U,1)),i);
    u(2,:,i) = U((0.5*size(U,1)+1):end,i);
end

% Use jth loading
j = loading;
c={'b','g','r','m','k','c'};
c{7} = [203 108 57] ./ 255;  % brown
if(disp)
    figure(1),clf;
end
for i=-3:3
    q{i+4} = avgq{1,1} + 0.5*i*sqrt(S(j,j))*u(:,:,j);
    if(Closed)
        q{i+4} = ProjectC(q{i+4});
    end
    X{i+4} = q_to_curve(q{i+4});
    X{i+4} = X{i+4} - repmat(mean(X{i+4},2),1,size(X{i+4},2));
    if(disp)
        hold on
        plot(X{i+4}(1,:)+0.23*i,X{i+4}(2,:),'Color',c{i+4},'LineWidth',1.5);
    end
    [n,N] = size(X{i+4});
    for m = 1:n
        z(m,:) = gradient(X{i+4}(m,:),1/(N-1));
    end
    len(i+4) = trapz(linspace(0,1,N),sqrt(sum(z.*z)));
end

if(disp)
    axis equal;
    %legend('-3','-2','-1','average','1','2','3')
    axis off;
end
