function [ts,vals,regions,numRegions,regionMaxSignificances] =  ...
            findRegionDensityTimeSeries(zVals,Q,xx,tSpacing,tWidth,...
                                        sigMap,sigma,xRange,minRegionSize)


    if nargin < 4 || isempty(tSpacing)
        tSpacing = 20;
    end
    
    if nargin < 5 || isempty(tWidth)
        tWidth = 50;
    end
    
    if nargin  < 7 || isempty(sigma)
        sigma = 1.5;
    end
    
    if nargin  < 8 || isempty(xRange) || length(xRange) ~= 2
        xRange = [-105 105];
    end
    
    if nargin  < 9 || isempty(minRegionSize)
        minRegionSize = 5;
    end
    
 
    if ~iscell(Q)
        
        s = size(Q);
        
        regions = bwlabel(Q);
        regionValues = setdiff(unique(regions),0);
        numRegions = length(regionValues);
        
        regionMaxSignificances = zeros(numRegions,1);
        for i=1:numRegions
            regionMaxSignificances(i) = min(sigMap(regions == regionValues(i)));
        end
        
        [~,sortIdx] = sort(regionMaxSignificances);
        regionMaxSignificances = regionMaxSignificances(sortIdx);
        temp = uint8(size(regions));
        for i=1:length(sortIdx)
            temp(regions == regionValues(sortIdx(i))) = i;
        end
        regionValues = setdiff(unique(temp),0);
        regions = temp;
        
        
        regionMaps = false(s(1),s(2),numRegions);
        useRegions = true(numRegions,1);
        for i=1:numRegions
            regionMaps(:,:,i) = regions == regionValues(i);
            if sum(sum(regionMaps(:,:,i))) < minRegionSize
                useRegions(i) = false;
                regions(regions == regionValues(i)) = 0;
            end
        end
        regionMaps = regionMaps(:,:,useRegions);
        numRegions = sum(useRegions);
        a = setdiff(unique(regions),0);
        for i=1:length(a)
            regions(regions == a(i)) = i;
        end
        
        numPoints = s(1);
        N = length(zVals);
        L = max(returnCellLengths(zVals))/2;
        ts = tSpacing:tSpacing:L;
        tt = -tWidth:tWidth;
        
        vals = zeros(12,length(ts),numRegions);
        parfor i=1:N
            fprintf('\t Calculating for Data Set #%2i out of %2i\n',i,N);
            temp = zeros(length(ts),numRegions);
            for j=1:length(ts);
                q = ts(j) + tt;
                q = q(q > 0 & q<= length(zVals{i}));
                [~,dd] = findPointDensity(zVals{i}(q,:),sigma,numPoints,xRange);
                for k=1:numRegions
                    temp(j,k) = sum(dd(regionMaps(:,:,k)))*(xx(2)-xx(1))^2;
                end
            end
            vals(i,:,:) = temp;
        end
        
    else
    
        s = size(Q{1});
        M = length(Q);
        
        numPoints = s(1);
        N = length(zVals);
        L = max(returnCellLengths(zVals))/2;
        ts = tSpacing:tSpacing:L;
        tt = -tWidth:tWidth;
        
        
        regions = cell(M,1);
        regionValues = cell(M,1);
        regionMaps = cell(M,1);
        numRegions = zeros(M,1);
        regionMaxSignificances = cell(M,1);
        
        for i=1:M
            
            regions{i} = bwlabel(Q{i});
            regionValues{i} = setdiff(unique(regions{i}),0);
            numRegions(i) = length(regionValues{i});
            
            regionMaxSignificances{i} = zeros(numRegions(i),1);
            for j=1:numRegions(i)
                regionMaxSignificances{i}(j) = max(sigMap(regions{i} == regionValues{i}(j)));
            end
            
            [~,sortIdx] = sort(regionMaxSignificances{i});
            regionMaxSignificances{i} = regionMaxSignificances{i}(sortIdx);
            
            temp = uint8(zeros(size(regions{i})));
            for j=1:length(sortIdx)
                temp(regions{i} == regionValues{i}(sortIdx(j))) = j;
            end
            regionValues{i} = setdiff(unique(temp),0);
            regions{i} = temp;
            
            
            regionMaps{i} = false(s(1),s(2),numRegions(i));
            
            useRegions = true(numRegions(i),1);
            for j=1:numRegions(i)
                regionMaps{i}(:,:,j) = regions{i} == regionValues{i}(j);
                if sum(sum(regionMaps{i}(:,:,j))) < minRegionSize
                    useRegions(j) = false;
                    regions{i}(regions{i} == regionValues{i}(j)) = 0;
                end
            end
            
            regionMaps{i} = regionMaps{i}(:,:,useRegions);
            numRegions(i) = sum(useRegions);
            a = setdiff(unique(regions{i}),0);
            for j=1:length(a)
                regions{i}(regions{i} == a(j)) = j;
            end
            
            
        end
          
        
        vals = zeros(12,length(ts),sum(numRegions));
        cSums = [0; cumsum(numRegions)];
        dx = xx(2) - xx(1);
        parfor i=1:N
            fprintf('\t Calculating for Data Set #%2i out of %2i\n',i,N);
            temp = zeros(length(ts),sum(numRegions));
            for j=1:length(ts);
                q = ts(j) + tt;
                q = q(q > 0 & q<= length(zVals{i}));
                [~,dd] = findPointDensity(zVals{i}(q,:),sigma,numPoints,xRange);
                for m=1:M
                    for k=1:numRegions(m)
                        r = cSums(m) + k;
                        temp(j,r) = sum(dd(regionMaps{m}(:,:,k)))*dx^2;
                    end
                end
            end
            vals(i,:,:) = temp;
        end
    
    end
    
    vals(isnan(vals) | isinf(vals)) = 0;
        
    
    
    
    
    