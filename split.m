% Splits curves into pieces using landmarks
% input: X - cell array of (2 x n) curves
% return: X - cell array of reorganized curves
% return: m - cell array of curve pieces

function [X,m,s] = split(X)

[X,s] = cursor(X);
numCurves = size(s,2);
numLm = length(s{1,1});

for i = 2:numCurves
    if length(s{1,i}) ~= length(s{1,i-1})
        error('Number of landmarks must be the same for each curve.');
    end
end

m = cell(numCurves,numLm+1);

for j = 1:numCurves
    p = X{1,j};
    if numLm == 0
        m{j,1} = p;
%         ql{j,1}=qp;
    else
        beg = 1;
        for k = 1:numLm
            m{j,k} = p(:,beg:s{1,j}(k));
            beg = s{1,j}(k);
        end

        m{j,numLm+1} = p(:,s{1,j}(numLm):end);
    end
end