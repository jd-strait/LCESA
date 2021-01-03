% Landmark selection for closed curves
% input: X - cell array of (2 x n) curves
% return: X - cell array of reorganized curves
% return: s - cell array of landmark indices

function [X,s] = cursor(X)
numCurves = size(X,2);
close all;
s = cell(1,numCurves);

% numCurves typically 2 if we are comparing 2 closed curves
% This loop plots each curve, asks to select a starting point and landmarks
for i = 1:numCurves
    
    p=X{1,i};
    figure(1),clf,hold on;
    plot(p(1,:),p(2,:),'LineWidth',3);
    
    h = datacursormode;
    set(h, 'DisplayStyle','datatip','SnapToDataVertex','on');
    disp('Select starting point');
    axis equal;
    v=axis;
    axis(v);
    pause;
    c = getCursorInfo(h);
    begInd = c(1).DataIndex;
    
    % This rearranges points of curve to start at the starting point
    if begInd ~= 1
        p = [p(:,begInd:(end-1)),p(:,1:begInd)];
    end
    
    % Plots first point, with starting point marked as red circle
    figure(1),clf,hold on;
    plot(p(1,:),p(2,:),'LineWidth',3);
    plot(p(1,1),p(2,1),'ro','LineWidth',3);
    h = datacursormode;
    set(h, 'DisplayStyle','datatip','SnapToDataVertex','on');
    disp('Select landmarks');
    axis(v);
    pause;
    c = getCursorInfo(h);
    
    figure(1),clf,hold on;
    plot(p(1,:),p(2,:),'LineWidth',3);
    plot(p(1,1),p(2,1),'ro','Linewidth',3);
    axis off;
    
    % numPoints is total number of points selected (starting + landmarks)
    numPoints = size(c,2);
    s{1,i} = zeros(numPoints,1);
    
    for j = 1:numPoints
        s{1,i}(j) = c(j).DataIndex;
        coord = c(j).Position;
        plot(coord(1),coord(2),'go','LineWidth',3);
    end
    axis(v);
    axis off;
    s{1,i} = sort(s{1,i});
    X{1,i} = p;
    pause;
end

close all;