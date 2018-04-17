function outputData = makeStimulusTriggeredPlots(ts,vals,LEDs,isControl,numSigma,timeBins,yBins,t1,t2,spacing)

    if nargin < 10 || isempty(spacing)
        spacing = 60;
    end
    t_min = -spacing/2 + .1;
    t_max = spacing/2 - .1;
    
    N = length(LEDs);
    W = zeros(N,length(ts));
    counts = zeros(N,1);
    for j=1:N
        CC = bwconncomp(LEDs{j});
        qq = returnFirstCellEntries(CC.PixelIdxList);
        for i=1:length(ts)
            r = argmin(abs(ts(i) - qq));
            W(j,i) = ts(i) - qq(r);
        end
        counts(j) = CC.NumObjects;
    end
    
    
    tVals = linspace(t_min,t_max,timeBins);  
    [~,~,bins] = histcounts(W./100,tVals);
    
    idx = find(~isControl);
    yy = linspace(0,1,yBins);
    Z = zeros(timeBins,yBins);
    for i=1:timeBins
        for k=1:length(idx)
            j = idx(k);
            q = vals(j,bins(j,:)==i);
            Z(i,:) = Z(i,:) + hist(q(:),yy);
        end
    end
    Z = bsxfun(@rdivide,Z,sum(Z,2));
    
    
    idx = find(isControl);
    yy = linspace(0,1,yBins);
    Z2 = zeros(timeBins,yBins);
    for i=1:timeBins
        for k=1:length(idx)
            j = idx(k);
            q = vals(j,bins(j,:)==i);
            Z2(i,:) = Z2(i,:) + hist(q(:),yy);
        end
    end
    Z2 = bsxfun(@rdivide,Z2,sum(Z2,2));
    
    
    m = sum(bsxfun(@times,yy,Z),2)./sum(yy);
    s = sqrt(sum(bsxfun(@minus,bsxfun(@times,yy,Z),m).^2,2)./sum(yy))/sqrt(sum(counts(~isControl))-1);
    
    m2 = sum(bsxfun(@times,yy,Z2),2)./sum(yy);
    s2 = sqrt(sum(bsxfun(@minus,bsxfun(@times,yy,Z2),m2).^2,2)./sum(yy))/sqrt(sum(counts(isControl))-1);
    
    X = [tVals(1:end-1) fliplr(tVals(1:end-1)) tVals(1)];
    Y = [m(1:end-1)-numSigma*s(1:end-1); flipud(m(1:end-1)+numSigma*s(1:end-1));m(1)-numSigma*s(1)];
    Y2 = [m2(1:end-1)-numSigma*s2(1:end-1); flipud(m2(1:end-1)+numSigma*s2(1:end-1));m2(1)-numSigma*s2(1)];
    
    plot(tVals,m,'r-','linewidth',2)
    hold on
    plot(tVals,m2,'b-','linewidth',2)
    fill(X,Y,'r','edgealpha',0,'facealpha',.5)   
    fill(X,Y2,'b','edgealpha',0,'facealpha',.5)
    
    set(gca,'fontsize',14,'fontweight','bold')
    xlabel('Time from Stimulus Onset (s)','fontsize',16)
    ylabel('Density in Region','fontsize',16)
    xlim([-30 30])
    q = ylim; 
    q(1) = 0;
    plot([t1 t1],[0 q(2)],'k--')
    plot([t2 t2],[0 q(2)],'k--')
    ylim(q)
    x = [t1 t2 t2 t1 t1];
    y = q([1 1 2 2 1]);
    %fill(x,y,[0 1 0],'facealpha',.2,'edgealpha',0);
    %legend({'Exp','Control'},'fontsize',14,'fontweight','bold')
    
    outputData.X = X;
    outputData.Y = Y;
    outputData.Y = Y2;
    outputData.m = m;
    outputData.s = s;
    outputData.m2 = m2;
    outputData.s2 = s2;
    outputData.ts = tVals;
    
    ylim([0 q(2)]);
    
    
    
    
    