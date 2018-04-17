function [ts,vals,LED_On_Off] = ...
    findChrimsonTimeSeries_densities_setTimes(zValues,LED,xx,...
                                on_window,off_window,sigma,sortOnOff)

                            
    if nargin < 4 || isempty(on_window)
        on_window = 0:1500;
    end
    
    if nargin < 5 || isempty(off_window)
        %should be between 0 and 6000
        off_window = 3000:4500;
    end
    
    if nargin < 6 || isempty(sigma)
        sigma = 1.5;
    end

    if nargin < 7 || isempty(sortOnOff)
        sortOnOff = true;
    end

       
    xRange = [min(xx) max(xx)];
    numPoints = length(xx);
    
    L = length(zValues(:,1)); 
    CC = bwconncomp(LED(1:L));
    zeroIdx = findFirstElementsInCells(CC.PixelIdxList);
    
    N = length(zeroIdx);
    onRegions = cell(N,1);
    offRegions = cell(N+1,1);
    for i=1:N
        onRegions{i} = zeroIdx(i) + on_window;
        onRegions{i} = onRegions{i}(onRegions{i} > 0 & onRegions{i} <= L);
        offRegions{i} = zeroIdx(i) - 6000 + off_window;
        offRegions{i} = offRegions{i}(offRegions{i} > 0 & offRegions{i} <= L);
    end
    
    offRegions{N+1} = zeroIdx(i) + off_window;
    offRegions{N+1} = offRegions{N+1}(offRegions{N+1} > 0 & offRegions{N+1} <= L);
    
    onRegions = onRegions(returnCellLengths(onRegions) > 0);
    offRegions = offRegions(returnCellLengths(offRegions) > 0);
    L_on = length(onRegions);
    L_off = length(offRegions);
    
    num = L_on + L_off;
    vals_on = zeros(numPoints,numPoints,L_on);
    vals_off = zeros(numPoints,numPoints,L_off);
    ts_on = zeros(L_on,1);
    ts_off = zeros(L_off,1);
    
    for i=1:length(onRegions)
        
        [~,d] = findPointDensity(zValues(onRegions{i},:),sigma,length(xx),xRange);
        vals_on(:,:,i) = d;
        ts_on(i) = round(median(onRegions{i}));
        
    end
    
    
    for i=1:length(offRegions)
        
        [~,d] = findPointDensity(zValues(offRegions{i},:),sigma,length(xx),xRange);
        vals_off(:,:,i) = d;
        ts_off(i) = round(median(offRegions{i}));
        
    end
    
    ts = [ts_on;ts_off];
    vals = zeros(numPoints,numPoints,num);
    vals(:,:,1:L_on) = vals_on;
    vals(:,:,L_on+1:end) = vals_off;
    LED_On_Off = [true(L_on,1);false(L_off,1)];
    
    if sortOnOff
        [~,sortIdx] = sort(ts);
        ts = ts(sortIdx);
        vals = vals(:,:,sortIdx);
        LED_On_Off = LED_On_Off(sortIdx);
    end
    
    
    
    
    
    
    
    