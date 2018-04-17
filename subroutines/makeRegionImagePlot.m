function outImage = makeRegionImagePlot(data,f,vecs,pixels,thetas,imageSize,baselineVals)

    d = length(data(1,:));
    L = length(f);
    numModes = d / L;
    numPixels = length(pixels);
   
    if nargin < 7 || isempty(baselineVals)
        baselineVals = zeros(numModes,1);
    end
        
    powerInMode = zeros(numModes,1);
    for i=1:numModes
        powerInMode(i) = sum(data((1:L) + (i-1)*L));
    end
    powerInMode(1:2) = 0;
    powerInMode = powerInMode ./ sum(powerInMode);
    
    outVals = zeros(numPixels,1);
    compVals = rectify(powerInMode-baselineVals);
    for i=1:numModes
        outVals = outVals + abs(vecs(:,i)).*compVals(i);
    end
    
    outRadonImage = zeros(imageSize);
    outRadonImage(pixels) = outVals;
    
    outImage = iradon(outRadonImage,thetas);