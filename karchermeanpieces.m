% Each row needs to be one of the curves with same number of pieces across
% all curves
function [nvbar,avg,totdist] = karchermeanpieces(X)
    numCurves = size(X,1);
    numPieces = size(X,2);
    N = 100; M = 300;
    maxIter = 200;
    eps = 0.05;
    nvbar(1)= 1;
    k = 2;
    totdist = 0;
    n = size(X{1,1},1);

for i=1:numCurves
    for j=1:numPieces
        X{i,j} = ReSampleCurve(X{i,j},N);
        X{i,j} = ReSampleCurve(X{i,j},N);
    end
end

for i=1:numCurves
    for j=1:numPieces-1
        Y{1,i}(:,(N*j - N - j + 2):((N-1)*j)) = X{i,j}(:,1:(N-1));
        Y{1,i}(:,(N*j - N - j + 2):((N-1)*j)) = X{i,j}(:,1:(N-1));
    end
    Y{1,i}(:,(N*numPieces - N - numPieces + 2):(((N-1)*numPieces))+1) = X{i,numPieces}(:,1:N);
    Y{1,i}(:,(N*numPieces - N - numPieces + 2):(((N-1)*numPieces))+1) = X{i,numPieces}(:,1:N);
end
    
for i = 1:numCurves
    % Resample and center curve
    Y{1,i} = ReSampleCurve(Y{1,i},300);
    cenmass(:,i) = mean(Y{1,i},2);
    Y{1,i} = Y{1,i} - repmat(cenmass(:,i),1,300);
    
    % Rescale curve
    for j = 1:n
        v(j,:) = gradient(Y{1,i}(j,:),1/(300-1));
    end
    len1(i) = trapz(linspace(0,1,300),sqrt(sum(v.*v)));
    Y{1,i} = Y{1,i}/len1(i);
    
    q{1,i} = curve_to_q(Y{1,i});
    q{1,i} = ProjectC(q{1,i});
    Y{1,i} = q_to_curve(q{1,i});
end

for i=1:numCurves
    for j=1:numPieces
        X{i,j} = X{i,j} - repmat(cenmass(:,i),1,size(X{i,j},2));
        X{i,j} = X{i,j}/len1(i);
    end
end
    
Nnew=N*numPieces - (numPieces-1);   
    
    % initialize on a circle
%     avgt = linspace(0,2*pi,N*numPieces - (numPieces-1));
%     avgcurve{1,1}(1,:) = cos(avgt);
%     avgcurve{1,1}(2,:) = sin(avgt);
%     avgcurve{1,1} = avgcurve{1,1}(:,end:-1:1);

    % Resample and center average
%     avgcurve{1,1} = ReSampleCurve(avgcurve{1,1},N);
%     avgcurve{1,1} = avgcurve{1,1} - repmat(mean(avgcurve{1,1},2),1,N);
        
    % Rescale average
%     len1 = sqrt(trapz(linspace(0,1,N),sum(gradient(avgcurve{1,1},1/(N)).^2,1)));
%     avgcurve{1,1} = avgcurve{1,1}/len1;
%     qc = curve_to_q(avgcurve{1,1});
%     avgcurve{1,1} = q_to_curve(qc);
%     [avgcurve,avgcurvel] = split(avgcurve);
%     for j=1:numPieces
%         avgcurvel{1,j}=ReSampleCurve(avgcurvel{1,j},N);
%         avgql{1,j}=curve_to_q2(avgcurvel{1,j});
%     end
%     
%     avgcurve{1,1}=avgcurvel{1,1};
%     for j=2:numPieces
%         avgcurve{1,1} = [avgcurve{1,1},avgcurvel{1,j}(:,2:N)];
%     end
%     avgq{1,1} = curve_to_q(avgcurve{1,1});

    len2 = repmat(0,numCurves,1);
    for i=1:numCurves
        for j=1:numPieces
            X{i,j} = ReSampleCurve(X{i,j},N);
            ql{i,j} = curve_to_q2(X{i,j});
            len2(i) = len2(i) + InnerProd_Q(ql{i,j},ql{i,j});
        end
        for j=1:numPieces
            ql{i,j} = ql{i,j}/sqrt(len2(i));
        end
    end
    
    for i=1:numCurves
        Y{1,i}=X{i,1};
        for j=2:numPieces
            Y{1,i} = [Y{1,i},X{i,j}(:,2:N)];
        end
        q{1,i} = curve_to_q(Y{1,i});
    end
    
    %initialize with one of the curves
    A = ceil(rand*size(Y,2));
    avgcurve{1,1}=Y{1,A};
    avgq{1,1}= q{1,A};
    for j=1:numPieces
        avgcurvel{1,j}=X{A,j};
        avgql{1,j}=ql{A,j};
%             avgcurvel{1,j} = ReSampleCurve(avgcurvel{1,j},N);
%             avgql{1,j} = curve_to_q2(avgcurvel{1,j});
    end
    stopcrit=0;
    
    while nvbar(k-1) > 0.03 && k <= maxIter && stopcrit<0.97
        
        if (nvbar(k-1)<0.05)
            eps = .02;
        end
        vbar=zeros(2,Nnew);
        totdist(k) = 0;
        for i = 1:numCurves
            if (trapz(linspace(0,1,Nnew),sum((avgq{1,1}-q{1,i}).^2))<.000001)
                
                v=zeros(2,Nnew);
                
            else
            
            % a = X(i,:);
            [~,qn,dist] = align2(avgq{1,1},avgql,Y{1,i},ql(i,:));
            totdist(k) = totdist(k) + (dist^2);
            v = (dist/sin(dist))*(qn - (cos(dist)*avgq{1,1}));
            
            end

            vbar = vbar + v/numCurves;
        end
        
%         for j = 1:numPieces
%             vsum = 0;
%             for i = 1:numCurves
%                 vsum = vsum + v{i,j};
%             end
%             vbar{j} = vsum / numCurves;
%         end
%         
%         for j=1:numPieces
%             if j==numPieces
%                 vbar2(:,(N*j - N - j + 2):(((N-1)*j))+1) = vbar{j}(:,1:N);
%             else
%                 vbar2(:,(N*j - N - j + 2):((N-1)*j)) = vbar{j}(:,1:(N-1));
%             end
%         end;
     
        nvbar(k) = sqrt(InnerProd_Q(vbar,vbar));
        avgq{1,1} = cos(eps*nvbar(k))*avgq{1,1} + sin(eps*nvbar(k))*(vbar/nvbar(k));
        avgq{1,1} = ProjectC(avgq{1,1});
        
        avgql{1,1}=avgq{1,1}(:,1:N);
        for j=2:numPieces
            avgql{1,j} = avgq{1,1}(:,(j-1)*N-(j-2):j*N-(j-1));
        end
        
        totdist
        pause
        
        k = k+1;
        
        if (k>2)
            stopcrit=(totdist(1)-totdist(k-1))/totdist(1);
        end
        
        avg = q_to_curve(avgq{1,1});
        
        figure(3),clf;
        plot(avg(1,:),avg(2,:));
        axis equal;
        hold on;
        pause(0.1)
        disp(totdist);
        
    end