function runData = findTimeDensities(zValueFiles,LEDfile,isControl,...
                            cameraNumbers,xx,sigma,on_window,off_window)


    if nargin < 6 || isempty(sigma)
        sigma = 1;
    end
    
    if nargin < 7 || isempty(on_window)
        on_window = 0:1500;
    end
    
    if nargin < 8 || isempty(off_window)
        off_window = 3000:4500;
    end

           
    
    N = length(zValueFiles);
    
    if ischar(zValueFiles{1})
        zVals = cell(N,1);
        inConvHulls = cell(N,1);
        fracInHull = zeros(N,1);
        fprintf(1,'\t Loading Embedding Values\n');
        for i=1:N
            
            zValues = [];zGuesses = [];
            load(zValueFiles{i},'zValues','inConvHull','zGuesses');
            a = zValues;
            a(~inConvHull,:) = zGuesses(~inConvHull,:);
            zVals{i} = a;
            inConvHulls{i} = inConvHull;
            fracInHull(i) = mean(inConvHull);
            
            clear zValues inConvHull zGuesses
        end
        zValues = zVals;
        clear zVals;
        
    else
        
        inConvHulls = [];
        fracInHull = [];
        zValues = zValueFiles;
        
    end

    
    fprintf(1,'\t Finding LED Values\n');
    
    if ischar(LEDfile)
   
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
                idx2 = min([t(idx(CC.PixelIdxList{k}(end))),181000]);
                LEDs{i}(idx1:idx2) = true;
            end
            
        end
        
    else
        
        LEDs = LEDfile;
        
    end
   
    
    
    clear t idx1 idx2 CC data
    
    runData.LEDs = LEDs;
    
    
    fprintf(1,'\t Calculating Embedding Time Series\n');
    ts = cell(N,1);
    vals = cell(N,1);
    LED_On_Offs = cell(N,1);
    for i=1:N
        [ts{i},vals{i},LED_On_Offs{i}] = ...
            findChrimsonTimeSeries_densities_setTimes(zValues{i},...
            LEDs{i},xx,on_window,off_window,sigma,false);
    end
    
    
    runData.ts = ts;
    runData.vals = vals;
    runData.LED_On_Offs = LED_On_Offs;
    runData.isControl = isControl;
    runData.cameraNumbers = cameraNumbers;
    
