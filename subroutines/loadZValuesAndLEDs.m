function [zVals,LEDs,spacing] = loadZValuesAndLEDs(files,led_file,cameraNumbers,maxNumCycles)

    N = length(files);
    if nargin < 4 || isempty(maxNumCycles)
        maxNumCycles = -1;
    end
    
    zVals = cell(N,1);
    for i=1:N
        %fprintf(1,'\t Loading File #%2i out of %2i\n',i,N);
        load(files{i},'zValues','inConvHull','zGuesses');
        zVals{i} = zValues;
        zVals{i}(~inConvHull,:) = zGuesses(~inConvHull,:);
        clear zValues inConvHull zGuesses
    end
    
    
    LEDs = cell(N,1);
    lengths = returnCellLengths(zVals) ./ 2;
    data = importdata(led_file);
    led_vals = data.data(:,13);
    CC = bwconncomp(led_vals);
    truncate = false;
    if maxNumCycles > 0 && CC.NumObjects > maxNumCycles
        for j=maxNumCycles+1:CC.NumObjects
            led_vals(CC.PixelIdxList{j}) = 0;
            truncate = true;
        end
        firstZeroIdx = CC.PixelIdxList{maxNumCycles+1}(1) - 1;
        t_ends = data.data(firstZeroIdx,1:12);
        t_end = round(median(t_ends(t_ends>0)));
    end
    
    
    spacings = cell(N,1);
    
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
            idx2 = min([t(idx(CC.PixelIdxList{k}(end))),181000]);
            LEDs{i}(idx1:idx2) = true;
        end
        
        if truncate && t_end < length(zVals{i}(:,1))
            zVals{i} = zVals{i}(1:t_end,:);
        end
        
        LEDs{i} = LEDs{i}(1:length(zVals{i}(:,1)));
        
        CC = bwconncomp(LEDs{i});
        spacings{i} = diff(returnFirstCellEntries(CC.PixelIdxList));
        
    end
    
    spacing = median(cell2mat(spacings'));
    
    
    
    
    
    
    