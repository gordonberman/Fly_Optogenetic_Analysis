function entropies = findEntropyCurves(zValueFile,LEDfile,cameraNumbers,frames,...
                                        timeWidth,numBins,axisLimits,loadFiles)


                                    
    if nargin < 4 || isempty(frames)
        frames = -2250:3750;
    end
                                    
                                    
    if nargin < 5 || isempty(timeWidth)
        timeWidth = 20;
    end
    
    
    if nargin < 6 || isempty(numBins)
        numBins = 51;
    else
        numBins = round(numBins);
        if mod(numBins,2) ~= 1
            numBins = numBins + 1;
        end
    end
    
    
    if nargin < 7 || isempty(axisLimits)
        axisLimits = [-105 105 -105 105];
    else
        if length(axisLimits) == 2
            axisLimits = [axisLimits axisLimits];
        end
    end
    
    if nargin < 8 || isempty(loadFiles)
        loadFiles = true;
    end
    
    if loadFiles
        
        if ~iscell(zValueFile)
            
            load(zValueFile,'zValues','inConvHull','zGuesses')
            zValues(~inConvHull,:) = zGuesses(~inConvHull,:);
            
            clear inConvHull zGuesses
            
            lengths = length(zValues(:,1));
            data = importdata(LEDfile);
            
            led_vals = data.data(:,13);
            LEDs = false(lengths,1);
            t = data.data(:,cameraNumbers);
            z = find(t > 0);
            q = argmin(t(z));
            q = z(min(q));
            idx = find(t > 0);
            idx = idx(idx >= q);
            CC = bwconncomp(led_vals(idx));
            for k=1:CC.NumObjects
                idx1 = max([t(idx(CC.PixelIdxList{k}(1))),1]);
                idx2 = min([t(idx(CC.PixelIdxList{k}(end))),181000]);
                LEDs(idx1:idx2) = true;
            end
            
        else
            
            N = length(zValueFile);
            
            zVals = cell(N,1);
            fprintf(1,'\t\t Loading Embedding Values\n');
            for i=1:N
                load(zValueFile{i},'zValues','inConvHull','zGuesses');
                zValues(~inConvHull,:) = zGuesses(~inConvHull,:);
                zVals{i} = zValues;
                clear zValues inConvHull zGuesses
            end
            zValues = zVals;
            clear zVals;
            
            LEDs = cell(N,1);
            lengths = returnCellLengths(zValues) ./ 2;
            data = importdata(LEDfile);
            led_vals = data.data(:,13);
            
            for i=1:N
                LEDs{i} = false(lengths(i),1);
                t = data.data(:,cameraNumbers(i));
                z = find(t > 0);
                q = argmin(t(z));
                q = z(min(q));
                idx = find(t > 0);
                idx = idx(idx >= q);
                CC = bwconncomp(led_vals(idx));
                for k=1:CC.NumObjects
                    idx1 = max([t(idx(CC.PixelIdxList{k}(1))),1]);
                    idx2 = min([t(idx(CC.PixelIdxList{k}(end))),lengths(i)]);
                    LEDs{i}(idx1:idx2) = true;
                end
            end
            
            for i=1:length(LEDs)
                LEDs{i} = LEDs{i}(1:length(zValues{i}(:,1)));
            end
            
            
            zValues = cell2mat(zValues);
            LEDs = cell2mat(LEDs);
            
        end
        
    else
        
        for i=1:length(LEDfile)
            LEDfile{i} = LEDfile{i}(1:length(zValueFile{i}(:,1)));
        end
        zValues = cell2mat(zValueFile);
        LEDs = cell2mat(LEDfile);
        
    end
    
    
    CC = bwconncomp(LEDs);
    zeroIdx = returnFirstCellEntries(CC.PixelIdxList);
    
    N = length(frames);
    L = length(zeroIdx);
    xxs = cell(1,2);
    xxs{1} = linspace(axisLimits(1),axisLimits(2),numBins);
    xxs{2} = linspace(axisLimits(3),axisLimits(4),numBins);
    dx2 = (xxs{1}(2) - xxs{1}(1))*(xxs{2}(2) - xxs{2}(1));
    
    entropies = zeros(N,1);
    for ww=1:N
 
        i = frames(ww);

        currentIdx = int32(zeros(L*(2*timeWidth+1),1));
        for j=1:(2*timeWidth+1)
            currentIdx(L*(j-1) + (1:L)) = zeroIdx + i - timeWidth + j;
        end
                
        a = currentIdx >= 1 & currentIdx <= length(LEDs);
        currentIdx = currentIdx(a);

        Z = hist3(zValues(currentIdx,:),xxs);
        Z = Z./(sum(Z(:)).*dx2);
        
        Z2 = Z ./ sum(Z(:));
        Hvals = -Z2.*log2(Z2);
        entropies(ww) = sum(Hvals(~isinf(Hvals) & ~isnan(Hvals)));
        
    end
    
    
    
    