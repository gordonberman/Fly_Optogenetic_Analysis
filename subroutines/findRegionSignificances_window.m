function runData = findRegionSignificances_window(...
    zValueFiles,LEDfile,isControl,cameraNumbers,xx,sigma,on_window,...
    off_window,densityThreshold,sigAlpha)


    if nargin < 6 || isempty(sigma)
        sigma = 1;
    end
    
    if nargin < 7 || isempty(on_window)
        on_window = 0:1500;
    end
    
    if nargin < 8 || isempty(off_window)
        off_window = 3000:4500;
    end

    if nargin < 9 || isempty(densityThreshold)
        densityThreshold = 1e-7;
    end
    
    if nargin < 10 || isempty(sigAlpha)
        sigAlpha = .05;
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
    
    [~,density] = findPointDensity(cell2mat(zValues),...
        sigma,length(xx),[xx(1) xx(end)]);
    
    
    runData.zValues = zValues;
    runData.inConvHulls = inConvHulls;
    runData.fracInHull = fracInHull;
    runData.zValueFiles = zValueFiles;
    runData.LEDfile = LEDfile;
    runData.isControl = isControl;
    runData.sigma = sigma;
    runData.cameraNumbers = cameraNumbers;
    runData.xx = xx;
    runData.density = density;
    runData.densityThreshold = densityThreshold;
    runData.sigAlpha = sigAlpha;
    runData.numIndividuals = length(runData.isControl);
    runData.numControl = sum(runData.isControl);
    runData.numExp = sum(~runData.isControl);
    
    
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
    
    
    if ~isempty(inConvHulls)
        runData.fracInHull_LED_On = zeros(N,1);
        for i=1:N
            a = length(inConvHulls{i});
            runData.fracInHull_LED_On(i) = mean(inConvHulls{i}(LEDs{i}(1:a)));
        end
    else
        runData.fracInHull_LED_On = [];
    end
    
    
    clear t idx1 idx2 CC data
    
    runData.LEDs = LEDs;
    
    
    fprintf(1,'\t Calculating Embedding Time Series\n');
    ts = cell(N,1);
    vals = cell(N,1);
    LED_On_Offs = cell(N,1);
    L_ons = zeros(N,1);
    L_offs = zeros(N,1);
    for i=1:N
        [ts{i},vals{i},LED_On_Offs{i}] = ...
            findChrimsonTimeSeries_densities_setTimes(zValues{i},...
            LEDs{i},xx,on_window,off_window,sigma);
        L_ons(i) = sum(LED_On_Offs{i});
        L_offs(i) = sum(~LED_On_Offs{i});
    end
    
    s = [length(xx) length(xx)];
    runData.s = s;
    runData.vals = vals;
    runData.ts = ts;
    runData.LED_On_Offs = LED_On_Offs;
    runData.vals = vals;
    L_on = max(L_ons);
    L_off = max(L_offs);
    runData.L_on = L_on;
    runData.L_off = L_off;
    
    
    fprintf(1,'\t Finding On/Off Values\n');
    chis = nan(L_on,N,s(1),s(2));
    
    for i=1:N
        
        offVals = nan(L_offs(i),s(1),s(2));
        onVals = nan(L_ons(i),s(1),s(2));
        
        for j=1:s(1)
            for k=1:s(2)
                a = squeeze(vals{i}(j,k,~LED_On_Offs{i}));
                offVals(1:length(a),j,k) = a;
                b = squeeze(vals{i}(j,k,LED_On_Offs{i}));
                onVals(1:length(b),j,k) = b;
            end
        end
        

        
        if L_offs(i) > L_ons(i)
            chis(1:L_ons(i),i,:,:) = ...
                onVals - .5*(offVals(1:end-1,:,:) + offVals(2:end,:,:));
        else
            if L_offs(i) == L_ons(i);        
                chis(1:L_ons(i),i,:,:) = onVals - offVals;
            else
                chis(1:L_offs(i),i,:,:) = onVals(1:L_offs(i),:,:) - offVals;
            end
        end
        
    end
    
    
    runData.chis = chis;
    
    all_chis_exp = reshape(chis(:,~isControl,:,:),[L_on*sum(~isControl) s(1) s(2)]);
    all_chis_control = reshape(chis(:,isControl,:,:),[L_on*sum(isControl) s(1) s(2)]);
    
    runData.median_exp = zeros(s(1:2));
    runData.median_control = zeros(s(1:2));
    
    fprintf(1,'\t Finding Sign Rank Test Values (On/Off)\n');
    rankSums_exp = zeros(s);
    for i=1:s(1)
        for j=1:s(2)
            if density(i,j) > densityThreshold
                a = squeeze(all_chis_exp(:,i,j));
                rankSums_exp(i,j) = signrank(a(~isnan(a)));
                runData.median_exp(i,j) = median(a(~isnan(a)));
                b = squeeze(all_chis_control(:,i,j));
                runData.median_control(i,j) = median(b(~isnan(b)));
            else
                rankSums_exp(i,j) = 1;
            end
        end
    end
    runData.rankSums_exp = rankSums_exp;
    
    if sum(isControl) > 0
        fprintf(1,'\t Finding Rank Sum Test Values (Experimental/Control)\n');
        rankSums = zeros(s);
        for i=1:s(1)
            for j=1:s(2)
                if density(i,j) > densityThreshold
                    a = squeeze(all_chis_exp(:,i,j));
                    b = squeeze(all_chis_control(:,i,j));
                    rankSums(i,j) = ranksum(a(~isnan(a)),b(~isnan(b)));
                else
                    rankSums(i,j) = 1;
                end
            end
        end
        
        runData.rankSums = rankSums;
    else
        runData.rankSums = zeros(s);
    end
    
    
    Hvals = -density.*log2(density);
    runData.entropy = sum(Hvals(~isnan(Hvals) & ~isinf(Hvals))) * (xx(2) - xx(1))^2;
    runData.numComparisons = ceil(2.^runData.entropy);
    
    runData.alpha = 1 - (1-sigAlpha)^(1/runData.numComparisons);
    
    runData.mask = runData.rankSums < runData.alpha & runData.rankSums_exp < runData.alpha;
    
    
    runData.outputMap = double(runData.mask).*runData.median_exp;
    

    
    
    
    
    


