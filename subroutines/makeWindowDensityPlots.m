function densities = makeWindowDensityPlots(windows,zValueCell,LEDcell,isControl,xx,sigma,dt)


    zValues = cell2mat(zValueCell(~isControl));
    LEDs = cell2mat(LEDcell(~isControl));

    CC = bwconncomp(LEDs);
    zeroIdx = returnFirstCellEntries(CC.PixelIdxList);
    
    L = length(zeroIdx);
    numPoints = length(xx);
    xRange = [xx(1) xx(end)];
    M = length(windows(:,1));
    
    densities = zeros(numPoints,numPoints,M);
    for ww=1:M
 
        minT = floor(min(windows(ww,:))/dt);
        maxT = ceil(max(windows(ww,:))/dt);
        timeWidth = maxT - minT + 1;      

        currentIdx = int32(zeros(L*timeWidth,1));
        for j=1:timeWidth
            currentIdx(L*(j-1) + (1:L)) = zeroIdx + minT + j - 1;
        end
                
        a = currentIdx >= 1 & currentIdx <= length(LEDs);
        currentIdx = currentIdx(a);

        [~,densities(:,:,ww)] = findPointDensity(zValues(currentIdx,:),sigma,numPoints,xRange);
        
    end