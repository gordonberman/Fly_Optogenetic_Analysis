function outputData = returnLEDValues(outputStats,LED_OFF_start,LED_OFF_end,...
                    LED_ON_start,LED_ON_end,dilateSize,numShuffles,...
                    maxNumPartitions,partitionReplicates)

    if nargin < 2 || isempty(LED_OFF_start)
        LED_OFF_start = -1.5;
    end
    
    if nargin < 3 || isempty(LED_OFF_end)
        LED_OFF_end = -.5;
    end
    
    if nargin < 4 || isempty(LED_ON_start)
        LED_ON_start = 0;
    end
    
    if nargin < 5 || isempty(LED_ON_end)
        LED_ON_end = 1;
    end
    
    if nargin < 6 || isempty(dilateSize)
        dilateSize = 5;
    end
    
    if nargin < 7 || isempty(numShuffles)
        numShuffles = 1000;
    end
    
    if nargin < 8 || isempty(maxNumPartitions)
        maxNumPartitions = 6;
    end
    
    if nargin < 9 || isempty(partitionReplicates)
        partitionReplicates = 24;
    end
    
    
    
    sigma = outputStats.parameters.sigma;
    numRegions = outputStats.numRegions(1);
    regionMap = outputStats.regions{1};
    parameters = outputStats.parameters;
    sessionName = outputStats.sessionName;
    maxNumCycles = parameters.maxNumCycles;
    outputData.originalRegionMap = regionMap;
    outputData.region_boundary = outputStats.B;
    outputData.sessionName = sessionName;
    
    newRegionMap = uint8(zeros(size(regionMap)));
    for i=1:numRegions
        A = regionMap == i;
        A = imdilate(A,strel('disk',dilateSize));
        newRegionMap(A) = i;
    end
    regionMap = newRegionMap;
    
    
    xRange = parameters.xRange;
    numPoints = parameters.numPoints;
    xx = linspace(xRange(1),xRange(2),numPoints);
    dx = xx(2) - xx(1);
    outputData.dx = dx;
    outputData.xx = xx;
    dt = 1/parameters.Fs;
    
    outputData.sigma = sigma;
    outputData.numRegions = numRegions;
    outputData.regionMap = regionMap;
    outputData.maxNumPartitions = maxNumPartitions;
    outputData.partitionReplicates = partitionReplicates;
    outputData.numShuffles = numShuffles;
    outputData.dilateSize = dilateSize;
    outputData.LED_OFF_start = LED_OFF_start;
    outputData.LED_OFF_end = LED_OFF_end;
    outputData.LED_ON_start = LED_ON_start;
    outputData.LED_ON_end = LED_ON_end;

    
    
    files = findImagesInFolder(parameters.embedding_path,'.mat');
    led_files = findImagesInFolder(parameters.led_file_path,'.txt');
    fileData = findFileStrainNamesAndCameras(files,led_files);
        
    idx = find(strcmp(fileData.filming_session_names,sessionName));
    q = fileData.filming_session_numbers == idx;
    files = fileData.files(q);
    cameraNumbers = fileData.cameraNumbers(q);
    isControl = fileData.isControl(q);
    
    outputData.isControl = isControl;
    outputData.cameraNumbers = cameraNumbers;
    outputData.files = files;
    
    led_file = led_files{strcmp(fileData.led_filming_sessions,sessionName)};
    
    
    [zValues,LEDs,spacing] = loadZValuesAndLEDs(files,led_file,cameraNumbers,maxNumCycles);
    N = length(zValues);
    
    [~,density_exp] = findPointDensity(cell2mat(zValues(~isControl)),sigma,length(xx),[xx(1) xx(end)]);
    [~,density_control] = findPointDensity(cell2mat(zValues(isControl)),sigma,length(xx),[xx(1) xx(end)]);
    density_exp = density_exp / (sum(density_exp(:))*dx^2);
    density_control = density_control / (sum(density_control(:))*dx^2);
    outputData.density_exp = density_exp;
    outputData.density_control = density_control;
    
    startTime = LED_ON_start;
    endTime = LED_ON_end;
    startTime_off = LED_OFF_start;
    endTime_off = LED_OFF_end;
    on_window = (startTime/dt):(endTime/dt);
    off_window = round(startTime_off/dt):round(endTime_off/dt);
    outputData.spacing = spacing;
    
    density_data = findTimeDensities(...
        zValues,LEDs,isControl,cameraNumbers,xx,sigma,on_window,off_window);
    
    after_region_means = cell(N,1);
    after_region_modes = cell(N,1);
    after_region_sums = cell(N,1);
    after_region_assignments = cell(N,1);
    outputData.prior_densities = cell(N,1);
    for i=1:N
        L = floor(length(density_data.vals{i}(1,1,:))/2);
        after_region_means{i} = zeros(L,2);
        after_region_modes{i} = zeros(L,2);
        after_region_sums{i} = zeros(L,numRegions+1);
        outputData.prior_densities{i} = zeros(numPoints,numPoints,L);
        after_region_assignments{i} = zeros(L,1);
        for j=1:L
            outputData.prior_densities{i}(:,:,j) = density_data.vals{i}(:,:,j+L);
            a = density_data.vals{i}(:,:,j);
            maxVal = max(a(:));
            [ii,jj]  = find(a == maxVal,1,'first');
            after_region_modes{i}(j,:) = xx([ii jj]);
            meanX = sum(sum(bsxfun(@times,a,xx')))/sum(a(:));
            meanY = sum(sum(bsxfun(@times,a,xx)))/sum(a(:));
            after_region_means{i}(j,:) = [meanX meanY];
            for k=1:numRegions
                after_region_sums{i}(j,k+1) = sum(sum(a.*double(regionMap == k)))*dx*dx;
            end
            after_region_sums{i}(j,1) = 1 - sum(after_region_sums{i}(j,2:end));
            
            regionValue = regionMap(ii,jj);
            after_region_assignments{i}(j) = regionValue;
            
        end
    end
        
    outputData.after_region_means = after_region_means;
    outputData.after_region_modes = after_region_modes;
    outputData.after_region_sums = after_region_sums;
    outputData.after_region_assignments = after_region_assignments;
    
    
    %outputData.density_data = density_data;
    
    
    
    
    
    %calculations for experimental data sets
    idx = find(~isControl);
    assignments = cell2mat(after_region_assignments(idx));
    
    densities = outputData.prior_densities(idx);
    Q = length(densities);
    densities = cell2mat(reshape(densities,[1 1 Q]));
    sum_densities = sum(sum(densities,1),2)*dx^2;
    densities = bsxfun(@rdivide,densities,sum_densities);

    MI_data = calculateStimulusMutualInformation(densities,...
                            assignments,density_exp,numRegions,xx);
                        
    outputData.MI_exp = MI_data.MI;
    outputData.average_prior_regions_exp = MI_data.averageRegions;
    outputData.prob_after_region_assignments_exp = MI_data.region_probs;
    
    if numShuffles > 1
        
        q = zeros(numShuffles,1);
        parfor i=1:numShuffles
            r = assignments
            w = r(randperm(length(assignments)));
            MI_data = calculateStimulusMutualInformation(densities,...
                            w,density_exp,numRegions,xx);
            q(i) = MI_data.MI;
        end    
        outputData.MI_exp_shuffles = q;
 
    end
        
    if maxNumPartitions > 1    
        outputData.MI_estimate_exp = calculateExtrapolatedMutualInformations(...
            assignments,densities,numRegions,maxNumPartitions,partitionReplicates,xx);     
        %makeExtrapolatedMIPlot(outputData.MI_estimate_exp)     
    end
    
    
    
    if numRegions > 1
    
        q = assignments > 0;
        if sum(q) > numRegions + 1
            
            assignments = assignments(q);
            densities = densities(:,:,q);
            newDensity = mean(densities,3);
            MI_data_no_zero = calculateStimulusMutualInformation(densities,...
                assignments,newDensity,numRegions,xx);
            
            outputData.MI_exp_no_zero = MI_data_no_zero.MI;
            outputData.average_prior_regions_exp_no_zero = MI_data_no_zero.averageRegions;
            outputData.prob_after_region_assignments_exp_no_zero = MI_data_no_zero.region_probs;
            
            if numShuffles > 1
                
                q = zeros(numShuffles,1);
                parfor i=1:numShuffles
                    r = assignments
                    w = r(randperm(length(assignments)));
                    MI_data = calculateStimulusMutualInformation(densities,...
                        w,newDensity,numRegions,xx);
                    q(i) = MI_data.MI;
                end
                outputData.MI_exp_shuffles_no_zero = q;
                
            end
            
            if maxNumPartitions > 1
                outputData.MI_estimate_exp_no_zero = calculateExtrapolatedMutualInformations(...
                    assignments,densities,numRegions,maxNumPartitions,partitionReplicates,xx);
                %makeExtrapolatedMIPlot(outputData.MI_estimate_exp)
            end
            
        end
        
    end
    
    
    
    
    
    %calculations for control data sets
    idx = find(isControl);
    assignments = cell2mat(after_region_assignments(idx));
    densities = outputData.prior_densities(idx);
    Q = length(densities);
    densities = cell2mat(reshape(densities,[1 1 Q]));
    sum_densities = sum(sum(densities,1),2)*dx^2;
    densities = bsxfun(@rdivide,densities,sum_densities);
    
    
    MI_data = calculateStimulusMutualInformation(densities,...
                            assignments,density_control,numRegions,xx);
                        
    outputData.MI_control = MI_data.MI;
    outputData.average_prior_regions_control = MI_data.averageRegions;
    outputData.prob_after_region_assignments_control = MI_data.region_probs;
    
    if numShuffles > 1
        
        q = zeros(numShuffles,1);
        parfor i=1:numShuffles
            r = assignments
            w = r(randperm(length(assignments)));
            MI_data = calculateStimulusMutualInformation(densities,...
                            w,density_control,numRegions,xx);
            q(i) = MI_data.MI;
        end    
        outputData.MI_control_shuffles = q;
        
    end                        
 
    if maxNumPartitions > 1     
        outputData.MI_estimate_control = calculateExtrapolatedMutualInformations(...
            assignments,densities,numRegions,maxNumPartitions,partitionReplicates,xx);      
        %makeExtrapolatedMIPlot(outputData.MI_estimate_control)
    end
    
    if numRegions > 1
    
        q = assignments > 0;
        if sum(q) > numRegions + 1
            
            assignments = assignments(q);
            densities = densities(:,:,q);
            newDensity = mean(densities,3);
            MI_data_no_zero = calculateStimulusMutualInformation(densities,...
                assignments,newDensity,numRegions,xx);
            
            outputData.MI_control_no_zero = MI_data_no_zero.MI;
            outputData.average_prior_regions_control_no_zero = MI_data_no_zero.averageRegions;
            outputData.prob_after_region_assignments_control_no_zero = MI_data_no_zero.region_probs;

            if numShuffles > 1
                
                q = zeros(numShuffles,1);
                parfor i=1:numShuffles
                    r = assignments
                    w = r(randperm(length(assignments)));
                    MI_data = calculateStimulusMutualInformation(densities,...
                        w,newDensity,numRegions,xx);
                    q(i) = MI_data.MI;
                end
                outputData.MI_control_shuffles_no_zero = q;
                
            end
            
            if maxNumPartitions > 1
                outputData.MI_estimate_control_no_zero = calculateExtrapolatedMutualInformations(...
                    assignments,densities,numRegions,maxNumPartitions,partitionReplicates,xx);
                %makeExtrapolatedMIPlot(outputData.MI_estimate_exp)
            end
            
        end
        
    end
    
    
    