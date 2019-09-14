function outputStats = analyzeSession_onlyDensities(sessionName,startTime,endTime,parameters)
%Code to run chrimson behavioral space analysis on a single filming session
%Outputs:
%   outputStats -> struct containing analysis information (see
%                   outputVariables.txt for details)
%
%Inputs:
%   sessionName -> string of session name (i.e. 'ss02316_0206')
%   startTime -> beginning of test region (in s), input empty array to
%                bring-up GUI (default = [])
%   endTime -> end of test region (in s), input empty array to bring-up GUI
%                (default = [])
%   parameters -> struct containing non-default run parameters (see
%                setRunParameters for definitions, default = [])


    addpath('utilities');
    addpath('subroutines');
    load('saved_colormaps');
    
    if nargin < 2
        startTime = [];
    end
    
    if nargin < 3
        endTime = [];
    end
    
    if nargin < 4
        parameters = [];
    end
    
    parameters = setRunParameters(parameters);

    fprintf(1,'Loading Training Set Data\n');
    load(parameters.training_set_path,'B','watershedMap')
    BB = B;
    
    fprintf(1,'Loading Embedding Data\n');
    files = findImagesInFolder(parameters.embedding_path,'.mat');
    led_files = findImagesInFolder(parameters.led_file_path,'.txt');
    fileData = findFileStrainNamesAndCameras(files,led_files);
    
    outputStats.fileData = fileData;
    outputStats.parameters = parameters;
        
    idx = find(strcmp(fileData.filming_session_names,sessionName));
    q = fileData.filming_session_numbers == idx;
    files = fileData.files(q);
    cameraNumbers = fileData.cameraNumbers(q);
    isControl = fileData.isControl(q);
    if parameters.no_controls
        isControl(:) = false;
    end
    
    led_file = led_files{strcmp(fileData.led_filming_sessions,sessionName)};
    
    outputStats.files = files;
    outputStats.led_file = led_file;
    outputStats.cameraNumbers = cameraNumbers;
    outputStats.isControl = isControl;
    maxNumCycles = parameters.maxNumCycles;
    outputStats.B = B;
    outputStats.sessionName = sessionName;
    
    [zValues,LEDs,spacing] = loadZValuesAndLEDs(files,led_file,cameraNumbers,maxNumCycles);
    
    xRange = parameters.xRange;
    numPoints = parameters.numPoints;
    xx = linspace(xRange(1),xRange(2),numPoints);
    dt = 1/parameters.Fs;
    t1 = parameters.entropyStartTime;
    spacing = round(spacing*dt);
    t2 = spacing + t1;
    entropies_ts = t1:dt:t2;
    timeWidth = parameters.entropyTimeWidth;
    numBins = parameters.entropyNumBins;
    sigma = parameters.sigma;
    
    outputStats.xx = xx;
    outputStats.entropies_ts = entropies_ts;
    outputStats.LEDs = LEDs;
    timeBins = parameters.timeBins;
    outputStats.LED_spacing = spacing;
    
    if isempty(startTime) || isempty(endTime) || parameters.use_entropy_window || parameters.force_calculate_entropy
        
        fprintf(1,'Calculating Entropies\n');
        
        if sum(~isControl) > 0
            entropies_exp = findEntropyCurves(zValues(~isControl),LEDs(~isControl),...
                cameraNumbers(~isControl),entropies_ts./dt,timeWidth,numBins,xRange,false);
            outputStats.entropies_exp = entropies_exp;
        else
            outputStats.entropies_exp = [];
        end
        
        
        if sum(isControl) > 0
            entropies_control = findEntropyCurves(zValues(isControl),LEDs(isControl),...
                cameraNumbers(isControl),entropies_ts./dt,timeWidth,numBins,xRange,false);
            outputStats.entropies_control = entropies_control;
        else
            outputStats.entropies_control = [];
            entropies_control = [];
        end
        
        
        
        test = ~parameters.use_entropy_window && (isempty(startTime) || isempty(endTime));
        if test
            figure
            xVals = labelEntropyPeak(entropies_ts,entropies_exp,entropies_control);
        end     
        
        if isempty(startTime) && ~parameters.use_entropy_window
            startTime = round(xVals(1)/dt)*dt;
        end
        
        if isempty(endTime) && ~parameters.use_entropy_window
            endTime = round(xVals(2)/dt)*dt;
        end
        
    end
    
    if parameters.use_entropy_window
        
        q = -entropies_exp';
        
        z = zeros(size(entropies_ts));
        dt_entropy = entropies_ts(2) - entropies_ts(1);
        num_t_entropy = round(startTime/dt_entropy);
        z(1:num_t_entropy) = -1;
        
        window_length = startTime;
        conv = ifft(conj(fft(z)).*(fft(q)));
        conv(entropies_ts < 0 | entropies_ts > 15) = 0;
        idx = find(conv > 0);
        minIdx = idx(argmin(conv(idx)));
        startTime = entropies_ts(minIdx);
        endTime = startTime + window_length;
        
        if startTime < 0
            endTime = endTime - startTime;
            startTime = 0;
        end

        if endTime > 15
            startTime = startTime - (endTime - 15);
            endTime = 15;
        end
                
        outputStats.startTime = startTime;
        outputStats.endTime = endTime;
        
    else
        
        outputStats.startTime = startTime;
        outputStats.endTime = endTime;
        
    end
        

    totalTime = endTime - startTime;
    windows = zeros(4,2);
    windows(1,:) = [startTime endTime] - totalTime;
    windows(2,:) = [startTime endTime];
    windows(3,:) = [startTime endTime] + totalTime;
    windows(4,:) = parameters.off_window;
    
    numRegions = max(watershedMap(:));
    outputStats.windows = windows;
    a = true(size(isControl));
    a(1) = false;
    wD = makeWindowDensityPlots(windows,zValues,LEDs,a,xx,sigma,dt);
    s_w = size(wD);
    window_densities_on_experimental = zeros(s_w(1),s_w(2),sum(~isControl));
    window_densities_off_experimental = zeros(s_w(1),s_w(2),sum(~isControl));
    window_densities_on_control = zeros(s_w(1),s_w(2),sum(isControl));
    window_densities_off_control = zeros(s_w(1),s_w(2),sum(isControl));
    window_region_densities_on_experimental = zeros(sum(~isControl),numRegions);
    window_region_densities_off_experimental = zeros(sum(~isControl),numRegions);
    window_region_densities_on_control = zeros(sum(isControl),numRegions);
    window_region_densities_off_control = zeros(sum(isControl),numRegions);
    
    if isControl(1)
        window_densities_on_control(:,:,1) = wD(:,:,2);
        window_densities_off_control(:,:,1) = wD(:,:,4);
        control_count = 2;
        experiment_count = 1;
    else
        window_densities_on_experimental(:,:,1) = wD(:,:,2);
        window_densities_off_experimental(:,:,1) = wD(:,:,4);
        control_count = 1;
        experiment_count = 2;
    end
    
    for i=2:length(isControl)
        a = true(size(isControl));
        a(i) = false;
        wD = makeWindowDensityPlots(windows,zValues,LEDs,a,xx,sigma,dt);
        if isControl(i)
            window_densities_on_control(:,:,control_count) = wD(:,:,2);
            window_densities_off_control(:,:,control_count) = wD(:,:,4);
            control_count = control_count + 1;
        else
            window_densities_on_experimental(:,:,experiment_count) = wD(:,:,2);
            window_densities_off_experimental(:,:,experiment_count) = wD(:,:,4);
            experiment_count = experiment_count + 1;
        end
    end
    
    
    regionIdx = cell(numRegions,1);
    for i=1:numRegions
        regionIdx{i} = find(watershedMap == i);
    end
    
    for i=1:sum(isControl)
        a = window_densities_on_control(:,:,i);
        b = window_densities_off_control(:,:,i);
        for j=1:numRegions
            window_region_densities_on_control(i,j) = sum(a(regionIdx{j}))*(xx(2)-xx(1))^2;
            window_region_densities_off_control(i,j) = sum(b(regionIdx{j}))*(xx(2)-xx(1))^2;
        end
    end
    window_region_densities_on_control = ...
        bsxfun(@rdivide,window_region_densities_on_control,sum(window_region_densities_on_control,2));
    window_region_densities_off_control = ...
        bsxfun(@rdivide,window_region_densities_off_control,sum(window_region_densities_off_control,2));
    
    for i=1:sum(~isControl)
        a = window_densities_on_experimental(:,:,i);
        b = window_densities_off_experimental(:,:,i);
        for j=1:numRegions
            window_region_densities_on_experimental(i,j) = sum(a(regionIdx{j}))*(xx(2)-xx(1))^2;
            window_region_densities_off_experimental(i,j) = sum(b(regionIdx{j}))*(xx(2)-xx(1))^2;
        end
    end
     window_region_densities_on_experimental = ...
        bsxfun(@rdivide,window_region_densities_on_experimental,sum(window_region_densities_on_experimental,2));
    window_region_densities_off_experimental = ...
        bsxfun(@rdivide,window_region_densities_off_experimental,sum(window_region_densities_off_experimental,2));
    
    
    %     wD = ...
    %         makeWindowDensityPlots(windows,zValues,LEDs,isControl,xx,sigma,dt);
    %     s_w = size(wD);
    %
    %         if sum(isControl) > 0
    %         wD_control = ...
    %             makeWindowDensityPlots(windows,zValues,LEDs,~isControl,xx,sigma,dt);
    %         window_densities = zeros(s_w(1),s_w(2),s_w(3) + length(wD_control(1,1,:)));
    %         window_densities(:,:,1:s_w(3)) = wD;
    %         window_densities(:,:,s_w(3)+1:end) = wD_control;
    %     else
    %         window_densities = wD;
    %     end
    
    outputStats.window_densities_on_experimental = window_densities_on_experimental;
    outputStats.window_densities_off_experimental = window_densities_off_experimental;
    outputStats.window_densities_on_control = window_densities_on_control;
    outputStats.window_densities_off_control = window_densities_off_control;
    outputStats.window_region_densities_on_experimental = window_region_densities_on_experimental;
    outputStats.window_region_densities_off_experimental = window_region_densities_off_experimental;
    outputStats.window_region_densities_on_control = window_region_densities_on_control;
    outputStats.window_region_densities_off_control = window_region_densities_off_control;
    
    
    %     fprintf(1,'Finding Significant Regions\n');
    %     on_window = (startTime/dt):(endTime/dt);
    %     off_window = round(parameters.off_window(1)/dt):round(parameters.off_window(end)/dt);
    %     densityThreshold = parameters.densityThreshold;
    %     sigAlpha = parameters.sigAlpha;
    %
    %     significanceData = findRegionSignificances_window(zValues,LEDs,isControl,...
    %         cameraNumbers,xx,sigma,on_window,off_window,densityThreshold,sigAlpha);
    %
    %     outputStats.differenceMap = significanceData.outputMap;
    %     outputStats.alpha = significanceData.alpha;
    %     outputStats.rankSums = significanceData.rankSums;
    %     outputStats.mask = significanceData.mask;
    %     outputStats.median_exp = significanceData.median_exp;
    %     outputStats.median_control = significanceData.median_control;
    %     outputStats.rankSums_exp = significanceData.rankSums_exp;
    %     outputStats.numComparisons = significanceData.numComparisons;
    %     outputStats.sigAlpha = significanceData.sigAlpha;
    %     outputStats.densityThreshold = significanceData.densityThreshold;
    %     outputStats.density = significanceData.density;
    %     outputStats.numIndividuals = length(significanceData.isControl);
    %     outputStats.numControls = sum(significanceData.isControl);
    %     outputStats.numExp = sum(~significanceData.isControl);
    %     outputStats.entropy = significanceData.entropy;
    %
    %     positiveMask = outputStats.mask > 0 & outputStats.differenceMap > 0;
    %     positiveMask = imdilate(positiveMask,strel('disk',parameters.dilateRegionSize));
    %     positiveMask = imclose(positiveMask,strel('disk',parameters.closeRegionSize));
    %     positiveMask = imfill(positiveMask,'holes');
    %     outputStats.positiveMask = positiveMask;
    %
    %     negativeMask = outputStats.mask > 0 & outputStats.differenceMap < 0;
    %     negativeMask = imdilate(negativeMask,strel('disk',parameters.dilateRegionSize));
    %     negativeMask = imclose(negativeMask,strel('disk',parameters.closeRegionSize));
    %     negativeMask = imfill(negativeMask,'holes');
    %     outputStats.negativeMask = negativeMask;
    %
    %     if max(positiveMask(:)) > 0 || max(negativeMask(:)) > 0
    %
    %         masks = {positiveMask,negativeMask};
    %
    %         fprintf(1,'Finding Region Time Series\n');
    %         tSpacing = parameters.triggered_plot_time_spacing;
    %         tWidth = parameters.triggered_plot_time_width;
    %         minRegionSize = parameters.minRegionSize;
    %
    %         A = zeros(numPoints,numPoints,2);
    %         A(:,:,1) = significanceData.rankSums_exp;
    %         A(:,:,2) = significanceData.rankSums;
    %         sigMatrix = min(A,[],3);
    %
    %         [ts,vals,regions,numRegions,regionMaxSigs] = findRegionDensityTimeSeries(...
    %             zValues,masks,xx,tSpacing,tWidth,sigMatrix,sigma,xRange,minRegionSize);
    %
    %         outputStats.ts = ts;
    %         outputStats.vals = vals;
    %         outputStats.regions = regions;
    %         outputStats.regionMaxSigs = regionMaxSigs;
    %         outputStats.sigMatrix = sigMatrix;
    %
    %     else
    %
    %         outputStats.ts = [];
    %         outputStats.vals = [];
    %         outputStats.regions = [];
    %         numRegions = [];
    %         outputStats.sigMatrix = [];
    %         outputStats.regionMaxSigs = [];
    %
    %     end
    %
    %     outputStats.numRegions = numRegions;
    %
    %     if ~isempty(numRegions)
    %         cSums = [0; cumsum(numRegions)];
    %     end
    %
    %     numSigma = parameters.numSigma;
    %     yBins = parameters.yBins;
    %     dx = xx(2) - xx(1);
    %     totalNumRegions = sum(numRegions);
    %     if ~isempty(numRegions)
    %         outputStats.totalNumRegions = totalNumRegions;
    %         outputStats.numPositiveRegions = numRegions(1);
    %         outputStats.numNegativeRegions = numRegions(2);
    %     else
    %         outputStats.totalNumRegions = [0 0];
    %         outputStats.numPositiveRegions = 0;
    %         outputStats.numNegativeRegions = 0;
    %     end
    %     timeBinSize = outputStats.parameters.timeBinSize;
    %     bootstrapNum = outputStats.parameters.bootstrapNum;
    %     confBound = outputStats.parameters.confBound;
    %
    %
    %     if parameters.plotsOn && sum(numRegions) > 0
    %
    %         fprintf(1,'Making Plots\n');
    %         stimulusPlots = cell(totalNumRegions,1);
    %         count = 1;
    %         for k=1:length(numRegions)
    %
    %             for i=1:numRegions(k)
    %
    %                 figure
    %                 r = cSums(k) + i;
    %
    %                 subplot(1,2,2)
    %                 currentVals = vals(:,:,r);
    %                 stimulusPlots{count} = makeStimulusTriggeredPlots(ts,currentVals,...
    %                     LEDs,isControl,numSigma,timeBins,yBins,startTime,endTime,outputStats.LED_spacing);
    %                 count = count + 1;
    %
    %                 subplot(1,2,1)
    %
    %                 imagesc(xx,xx,regions{k} == i);
    %
    %                 hold on
    %                 axis equal tight off xy
    %
    %                 B = bwboundaries(regions{k} == i);
    %                 for j=1:length(B)
    %                     plot(xx(B{j}(:,2)),xx(B{j}(:,1)),'k-','linewidth',2);
    %                 end
    %
    %                 for j=1:length(BB)
    %                     plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',3)
    %                 end
    %
    %                 if k == 1
    %                     title(['+ Region #' num2str(i)],'fontsize',16,'fontweight','bold')
    %                     colormap(cc)
    %                 else
    %                     title(['- Region #' num2str(i)],'fontsize',16,'fontweight','bold')
    %                     colormap(cc3)
    %                 end
    %
    %                 freezeColors;
    %
    %
    %                 drawnow
    %
    %                 makeTimeSeriesPlot(ts,currentVals,LEDs,isControl,timeBinSize,bootstrapNum,confBound);
    %                 if k == 1
    %                     title(['+ Region #' num2str(i)],'fontsize',16,'fontweight','bold')
    %                     colormap(cc)
    %                 else
    %                     title(['- Region #' num2str(i)],'fontsize',16,'fontweight','bold')
    %                     colormap(cc3)
    %                 end
    %
    %                 drawnow
    %
    %             end
    %
    %         end
    %
    %
    %         figure
    %         a = round(windows.*100)/100;
    %
    %         for i=1:4
    %             if sum(isControl) > 0
    %                 subplot(2,4,i)
    %             else
    %                 subplot(1,4,i)
    %             end
    %             imagesc(xx,xx,window_densities(:,:,i))
    %             axis equal tight off xy
    %             hold on;
    %             for j=1:length(BB)
    %                 plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',2)
    %             end
    %             title([num2str(a(i,1)) 's to ' num2str(a(i,2)) 's'],'fontsize',14,'fontweight','bold')
    %             colormap(cc)
    %             caxis([0 4e-4])
    %         end
    %
    %         if sum(isControl) > 0
    %             for k=1:4
    %                 i = k+4;
    %                 subplot(2,4,i)
    %                 imagesc(xx,xx,window_densities(:,:,i))
    %                 axis equal tight off xy
    %                 hold on;
    %                 for j=1:length(BB)
    %                     plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',2)
    %                 end
    %                 title([num2str(a(k,1)) 's to ' num2str(a(k,2)) 's (control)'],'fontsize',14,'fontweight','bold')
    %                 colormap(cc)
    %                 caxis([0 4e-4])
    %             end
    %         end
    %
    %         colorVals = 'brkmgc';
    %         segs = {'-','--','-.'};
    %
    %         if sum(isControl) > 0 && ~isempty(numRegions)
    %             figure
    %             hold on
    %             hs = zeros(numRegions(1),1);
    %             for i=1:numRegions(1)
    %                 tVals = stimulusPlots{i}.ts;
    %                 m = stimulusPlots{i}.m;
    %                 s = stimulusPlots{i}.s;
    %                 m2 = stimulusPlots{i}.m2;
    %                 s2 = stimulusPlots{i}.s2;
    %
    %                 m = m - m2;
    %                 s = sqrt(s.^2 + s2.^2);
    %
    %                 X = [tVals(1:end-1) fliplr(tVals(1:end-1)) tVals(1)];
    %                 Y = [m(1:end-1)-s(1:end-1); flipud(m(1:end-1)+s(1:end-1));m(1)-s(1)];
    %
    %                 hs(i) = plot(tVals,m,[colorVals(mod(i-1,6)+1) segs{floor((i-1)/5)+1}],'linewidth',2);
    %                 fill(X,Y,colorVals(mod(i-1,6)+1),'edgealpha',0,'facealpha',.5)
    %
    %             end
    %             yvals = ylim;
    %             plot([0 0],yvals,'k--','linewidth',1)
    %             plot([endTime endTime],yvals,'k--','linewidth',1)
    %             plot(stimulusPlots{i}.ts,zeros(size(stimulusPlots{i}.ts)),'k--','linewidth',1)
    %             ylim(yvals);
    %             set(gca,'fontsize',14,'fontweight','bold')
    %             xlabel('Time From Stimulus (s)','fontsize',16,'fontweight','bold')
    %             ylabel('Change in Region Density Above Controls','fontsize',16,'fontweight','bold')
    %             a = repmat({'Region #'},numRegions(1),1);
    %             for i=1:length(a)
    %                 a{i} = [a{i} num2str(i)];
    %             end
    %             legend(hs,a,'fontsize',16,'fontweight','bold','location','northwest')
    %         end
    %
    %         if ~isempty(numRegions)
    %         figure
    %         hold on
    %         hs = zeros(numRegions(1),1);
    %         for i=1:numRegions(1)
    %             tVals = stimulusPlots{i}.ts;
    %             m = stimulusPlots{i}.m;
    %             s = stimulusPlots{i}.s;
    %
    %             X = [tVals(1:end-1) fliplr(tVals(1:end-1)) tVals(1)];
    %             Y = [m(1:end-1)-s(1:end-1); flipud(m(1:end-1)+s(1:end-1));m(1)-s(1)];
    %
    %             hs(i) = plot(tVals,m,[colorVals(mod(i-1,6)+1) segs{floor((i-1)/5)+1}],'linewidth',2);
    %             fill(X,Y,colorVals(mod(i-1,6)+1),'edgealpha',0,'facealpha',.5)
    %
    %         end
    %         yvals = ylim;
    %         plot([0 0],yvals,'k--','linewidth',1)
    %         plot([endTime endTime],yvals,'k--','linewidth',1)
    %         plot(stimulusPlots{i}.ts,zeros(size(stimulusPlots{i}.ts)),'k--','linewidth',1)
    %         ylim(yvals);
    %         set(gca,'fontsize',14,'fontweight','bold')
    %         xlabel('Time From Stimulus (s)','fontsize',16,'fontweight','bold')
    %         ylabel('Total Region Density','fontsize',16,'fontweight','bold')
    %         a = repmat({'Region #'},numRegions(1),1);
    %         for i=1:length(a)
    %             a{i} = [a{i} num2str(i)];
    %         end
    %         legend(hs,a,'fontsize',16,'fontweight','bold','location','northwest')
    %
    %
    %         figure
    %         hold on
    %         hs = zeros(numRegions(1),1);
    %         for i=1:numRegions(1)
    %
    %             regionSize = sum(sum(regions{1}==i))*dx^2;
    %
    %             tVals = stimulusPlots{i}.ts;
    %             m = stimulusPlots{i}.m ./ regionSize;
    %             s = stimulusPlots{i}.s ./ regionSize;
    %
    %             X = [tVals(1:end-1) fliplr(tVals(1:end-1)) tVals(1)];
    %             Y = [m(1:end-1)-s(1:end-1); flipud(m(1:end-1)+s(1:end-1));m(1)-s(1)];
    %
    %             hs(i) = plot(tVals,m,[colorVals(mod(i-1,6)+1) segs{floor((i-1)/5)+1}],'linewidth',2);
    %             fill(X,Y,colorVals(mod(i-1,6)+1),'edgealpha',0,'facealpha',.5)
    %
    %         end
    %         yvals = ylim;
    %         plot([0 0],yvals,'k--','linewidth',1)
    %         plot([endTime endTime],yvals,'k--','linewidth',1)
    %         plot(stimulusPlots{i}.ts,zeros(size(stimulusPlots{i}.ts)),'k--','linewidth',1)
    %         ylim(yvals);
    %         set(gca,'fontsize',14,'fontweight','bold')
    %         xlabel('Time From Stimulus (s)','fontsize',16,'fontweight','bold')
    %         ylabel('Average Region Density','fontsize',16,'fontweight','bold')
    %         a = repmat({'Region #'},numRegions(1),1);
    %         for i=1:length(a)
    %             a{i} = [a{i} num2str(i)];
    %         end
    %         legend(hs,a,'fontsize',16,'fontweight','bold','location','northwest')
    %
    %
    %         figure
    %
    %         toPlot = setdiff(unique(regions{1}),0);
    %         L = length(toPlot);
    %         N = ceil(sqrt(L));
    %         M = ceil(L/N);
    %         for i=1:L
    %             subplot(M,N,i)
    %
    %             Q = zeros(size(regions{1}));
    %             Q(regions{1} == i) = i;
    %
    %             imagesc(xx,xx,Q);
    %
    %             hold on
    %             axis equal tight off xy
    %
    %             B = bwboundaries(regions{1} == i);
    %             for j=1:length(B)
    %                 plot(xx(B{j}(:,2)),xx(B{j}(:,1)),'k-','linewidth',2);
    %             end
    %
    %             for j=1:length(BB)
    %                 plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',3)
    %             end
    %
    %             title(['Region #' num2str(i)],'fontsize',16,'fontweight','bold')
    %             colormap(cc6)
    %             caxis([0 8])
    %
    %         end
    %
    %         end
    %
    %
    %     end
    
    
    
    
