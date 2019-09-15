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
    load(parameters.training_set_path,'watershedMap')
    
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
    
    
        
    outputStats.window_densities_on_experimental = window_densities_on_experimental;
    outputStats.window_densities_off_experimental = window_densities_off_experimental;
    outputStats.window_densities_on_control = window_densities_on_control;
    outputStats.window_densities_off_control = window_densities_off_control;
    outputStats.window_region_densities_on_experimental = window_region_densities_on_experimental;
    outputStats.window_region_densities_off_experimental = window_region_densities_off_experimental;
    outputStats.window_region_densities_on_control = window_region_densities_on_control;
    outputStats.window_region_densities_off_control = window_region_densities_off_control;
    
    
    
    
    
    
    
