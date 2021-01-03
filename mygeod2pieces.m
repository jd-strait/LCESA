function [dist,X2n,q2n,X1n]=mygeod2pieces(Y1,Y2)

%input: two curves Y1 and Y2 as cells that are 2xn or 3xn (2D vs. 3D
%curve); cells represent the curve broken up into different pieces already
%because of pre-specified landmarks; each curve must have same number of
%pieces (cells)
%output: dist=distance, also will display when you run the program; X2n:
%optimally registered curve X2; q2n: same as X2n except in q-fun form; X1n:
%normalized curve X1;

% Find optimal reparameterization?
reparam = 1;
gd = 0;    % gradient-descent
dp = 1;    % piecewise dynamic programming

% What displays you want to see, set to 1 if you want to see figures
Disp_geodesic_between_the_curves = 1;
Disp_registration_between_curves = 0;
Disp_optimal_reparameterization = 0;

% Resample each piece of the curve to have N points
numPieces = size(Y1,2);
N = 100;

for j=1:numPieces
    Y1{j} = ReSampleCurve(Y1{j},N);
    Y2{j} = ReSampleCurve(Y2{j},N);
end

% X1 and X2 are the joined together pieces (the whole curve)
for j=1:numPieces-1
    X2(:,(N*j - N - j + 2):((N-1)*j)) = Y2{j}(:,1:(N-1));
    X1(:,(N*j - N - j + 2):((N-1)*j)) = Y1{j}(:,1:(N-1));
end

X2(:,(N*numPieces - N - numPieces + 2):(((N-1)*numPieces))+1) = Y2{numPieces}(:,1:N);
X1(:,(N*numPieces - N - numPieces + 2):(((N-1)*numPieces))+1) = Y1{numPieces}(:,1:N);

% Center curves, not really needed but good for display purposes
cenmass1 = mean(X1,2);
cenmass2 = mean(X2,2);
X1 = X1 - repmat(cenmass1,1,size(X1,2));
X2 = X2 - repmat(cenmass2,1,size(X2,2));

for j=1:numPieces
    Y1{j} = Y1{j} - repmat(cenmass1,1,size(Y1{j},2));
    Y2{j} = Y2{j} - repmat(cenmass2,1,size(Y2{j},2));
end

% Rescale the curves
len1=0;
len2=0;
for j=1:numPieces
    for i = 1:2
        v1(i,:) = gradient(Y1{j}(i,:),1/(N-1));
        v2(i,:) = gradient(Y2{j}(i,:),1/(N-1));
    end
    len1 = len1 + trapz(linspace(0,1,N),sqrt(sum(v1.*v1)));
    len2 = len2 + trapz(linspace(0,1,N),sqrt(sum(v2.*v2)));
end

X1 = X1/len1;
X2 = X2/len2;
for j=1:numPieces
    Y1{j}=Y1{j}/len1;
    Y2{j}=Y2{j}/len2;
end

% Landmark indices
for i=1:2
    I{i} = linspace(N,size(X1,2)-(N-1),numPieces-1)';
end

keyboard;

% Align pieces using align function
if reparam == 1 && gd == 1
    [X1nGD,q1nGD,X2nGD,q2nGD,gamGDI,energy1,energy2] = alignpiecesGD(X1,X2,I,length(X1));
    gamGD = [linspace(0,1,length(gamGDI));gamGDI];
end

if reparam == 1 && dp == 1
    [X1nDP,q1nDP,X2nDP,q2nDP,q1l,q2ln,gamDP] = alignpiecesDP(X1,X2,Y1,Y2);
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

dist = distDP;
q1n = q1nDP;
q2n = q2nDP;
X1n = X1nDP;

sprintf('The distance between the two curves is %0.3f',dist)
colorspecstr = ['b','m','k','r','g','b','k','m','g','r','g','b','k'];
if(Disp_geodesic_between_the_curves)
    if (size(X1n,1)==2)
        for t=1:7
            s = dist*(t-1)/6;
            PsiQ(:,:,t) = (sin(dist - s)*q1n + sin(s)*q2n)/sin(dist);
            %                 PsiQ(:,:,t) = ProjectC(PsiQ(:,:,t));
            PsiX(:,:,t) = q_to_curve(PsiQ(:,:,t));
        end
        figure(4); clf; axis equal; hold on;
        for t=1:7
            z = plot(0.17*t + PsiX(1,:,t), PsiX(2,:,t),'r-');
            set(z,'LineWidth',[2],'color',colorspecstr(t));
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
            z = plot3(0.45*t + PsiX(1,:,t), PsiX(2,:,t), PsiX(3,:,t),'r-');
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
        z = plot(X1n(1,:), X1n(2,:),'r');
        set(z,'LineWidth',1.5);
        axis off;
        hold on;
        z = plot(X2n(1,:), 0.15+X2n(2,:),'b-+');
        set(z,'LineWidth',[1.5]);
        N = size(X1n,2);
        for j=1:N/15
            i = j*15;
            plot([X1n(1,i) X2n(1,i)],[X1n(2,i) 0.15+X2n(2,i)], 'k');
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
    plot(gam(1,:),gam(2,:),'LineWidth',2)
    axis equal tight
end