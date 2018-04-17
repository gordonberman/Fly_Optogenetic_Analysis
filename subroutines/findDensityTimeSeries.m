function [ts,vals] =  findDensityTimeSeries(zVals,xx,LEDs,sigma,timeBins,tSpacing)



    if nargin  < 4 || isempty(sigma)
        sigma = 1.5;
    end
        
    if nargin  < 5 || isempty(timeBins)
        timeBins = 300;
    end
    
    if nargin  < 6 || isempty(tSpacing)
        tSpacing = 20;
    end
    
    
    numPoints = length(xx);
    xRange = [min(xx) max(xx)];
    N = min([length(zVals),length(LEDs)]);
    L = max(returnCellLengths(zVals))/2;
    ts = tSpacing:tSpacing:L;
    
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
    
    tVals = linspace(-29.9,29.9,timeBins);
    [~,~,bins] = histcounts(W./100,tVals);
    
    vals = zeros(numPoints,numPoints,timeBins);
    parfor i=1:timeBins
        
        temp = cell(N,1);
        for j=1:N
            temp{j} = zVals{j}(bins==i,:);
        end

        [~,vals(:,:,i)] = findPointDensity(cell2mat(temp),sigma,numPoints,xRange);
        
    end

    
    