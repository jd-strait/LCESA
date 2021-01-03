function [I,X1,X2,Y1,Y2]=mygeod2(X1,X2)

%input: two curves X1 and X2 as 2xn or 3xn (2D vs. 3D curve)
%output: dist=distance, also will display when you run the program; X2n:
%optimally registered curve X2; q2n: same as X2n except in q-fun form; X1n:
%normalized curve X1;

% Find optimal reparameterization?
reparam = 0;
gd = 0;    % gradient-descent
dp = 0;    % piecewise dynamic programming

% What displays you want to see, set to 1 if you want to see figures
Disp_geodesic_between_the_curves = 1;
Disp_registration_between_curves = 1;
Disp_optimal_reparameterization = 1;

% Resample the curves to have N points
N = 800;
X1 = ReSampleCurve(X1,N);
X2 = ReSampleCurve(X2,N);

% Center curves, not really needed but good for display purposes
X1 = X1 - repmat(mean(X1,2),1,size(X1,2));
X2 = X2 - repmat(mean(X2,2),1,size(X2,2));

% Rescale curves
n = size(X1,1);
v1 = zeros(n,N);
v2 = zeros(n,N);
for i = 1:n
    v1(i,:) = gradient(X1(i,:),1/(N-1));
    v2(i,:) = gradient(X2(i,:),1/(N-1));
end

len1 = sqrt(InnerProd_Q(v1,v1));
X1 = X1/len1;
len2 = sqrt(InnerProd_Q(v2,v2));
X2 = X2/len2;

% Plot the two curves and pause
% figure(1),clf;
% plot(X1(1,:),X1(2,:),'LineWidth',2)
% axis equal
% figure(2),clf;
% plot(X2(1,:),X2(2,:),'LineWidth',2)
% axis equal
% 
% pause

% Run landmark program
[X,Y,I] = split({X1,X2});
Y1 = Y(1,:)';
Y2 = Y(2,:)';
X1 = X{1};
X2 = X{2};

numPieces = size(Y1,1);

% Align pieces using align function
if reparam == 1 && gd == 1
    tic
    [X1nGD,q1nGD,X2nGD,q2nGD,gamGDI,energy1,energy2] = alignGD(X1,X2,I,N);
    gamGD = [linspace(0,1,N);gamGDI];
    toc
end

if reparam == 1 && dp == 1
    tic
    [X1nDP,q1nDP,X2nDP,q2nDP,q1l,q2ln,gamDP] = alignDP(X1,X2,Y1,Y2);
    toc
end

% Forming geodesic between the registered curves
if dp == 1
    numPieces = size(q1l,1);
    for j=1:numPieces
        innprod(j)=InnerProd_Q(q1l{j},q2ln{j});
    end
    distDP = acos(sum(innprod));
end

if gd == 1
    distGD = acos(InnerProd_Q(q1nGD,q2nGD));
end

keyboard;

% Change between GD and DP here before displaying geodesic
dist = distDP;
q1n = q1nDP;
q2n = q2nDP;
X1n = X1nDP;

if reparam == 0
    X1n = X1;
    X2n = X2;
    q1n = curve_to_q(X1n);
    q2n = curve_to_q(X2n);
    dist = acos(InnerProd_Q(q1n,q2n));
end

% Run this if switching between after running first time
clear PsiQ PsiX
    
%sprintf('The distance between the two curves is %0.3f',dist)
colorspecstr = ['b','m','k','r','g','b','k','m','g','r','g','b','k'];
if(Disp_geodesic_between_the_curves)
    if (size(X1n,1)==2)
        for t=1:7
            s = dist*(t-1)/6;
            PsiQ(:,:,t) = (sin(dist - s)*q1n + sin(s)*q2n)/sin(dist);
            PsiQ(:,:,t) = ProjectC(PsiQ(:,:,t));
            PsiX(:,:,t) = q_to_curve(PsiQ(:,:,t));
        end
        figure(4); clf; axis equal; hold on;
        for t=1:7
            z = plot(0.14*t + PsiX(1,:,t), PsiX(2,:,t),'r-');
            %scatter(0.19*t+PsiX(1,1:99:end,t),PsiX(2,1:99:end,t),'Filled',colorspecstr(t));
            set(z,'LineWidth',2,'color',colorspecstr(t));
        end
        axis off;
    else
        for t=1:7
            s = dist*(t-1)/6;
            PsiQ(:,:,t) = (sin(dist - s)*q1n + sin(s)*q2n)/sin(dist);
            PsiX(:,:,t) = q_to_curve(PsiQ(:,:,t));
            PsiX(:,end+1,t)=PsiX(:,1,t);
        end
        figure(4); clf; axis equal; hold on;
        for t=1:7
            z = plot3(0.25*t + PsiX(1,:,t), PsiX(2,:,t), PsiX(3,:,t),'r-');
            set(z,'LineWidth',[2],'color',[(t-1)/6 (t-1)/12 0]);
        end
        axis off;
    end
end

% Displaying the correspondence
if(Disp_registration_between_curves)
    X2n=PsiX(:,:,end);
    X1n=PsiX(:,:,1);
    if (size(X1,1)==2)
        figure(3); clf;
        hold on;
        z = plot(X1n(1,:), X1n(2,:),'r');
        set(z,'LineWidth',1.5);
        %scatter(X1n(1,1:99:end),X1n(2,1:99:end),'Filled','r');
        axis off;
         z = plot(X2n(1,:), 0.15+X2n(2,:),'b-+');
% z = plot(X2n(1,:), 0.15+X2n(2,:),'b');
        set(z,'LineWidth',1.5);
%         scatter(X2n(1,1:99:end),0.15+X2n(2,1:99:end),'Filled','b');
        N = size(X1n,2);
        for j=1:N/15
            i = j*15;
            plot([X1n(1,i) X2n(1,i)],[X1n(2,i) 0.15+X2n(2,i)], 'k');
        end
        a=1:99:length(X1n);
        for i=1:length(a)
            plot([X1n(1,a(i)) X2n(1,a(i))],[X1n(2,a(i)) 0.15+X2n(2,a(i))],'k');
        end
    else
        figure(3); clf;
        z = plot3(X1n(1,:), X1n(2,:), X1n(3,:),'r');
        set(z,'LineWidth',[1.5]);
        axis off;
        hold on;
        z = plot3(X2n(1,:), 0.15+X2n(2,:), X2n(3,:),'b-+');
        set(z,'LineWidth',1.5);
        N = size(X1n,2);
        for j=1:N/15
            i = j*15;
            plot3([X1n(1,i) X2n(1,i)],[X1n(2,i) 0.15+X2n(2,i)], [X1n(3,i) X2n(3,i)], 'k');
        end
    end
end

if(Disp_optimal_reparameterization)
    figure(100);clf;
    plot(gamDP(1,:),gamDP(2,:),'LineWidth',2)
    axis square
end