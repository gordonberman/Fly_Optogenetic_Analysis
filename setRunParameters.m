function parameters = setRunParameters(parameters) 
%Sets all run parameters for analyzeSession.  To set any parameter to a 
%non-default value, change the appropriate value in the struct you pass
%into this function (i.e. parameters.sigAlpha = .01).  All values not
%included in the inputted struct will be set to the default values listed
%below.

    if nargin < 1
        parameters = [];
    end
    
    %path to *_embedding.mat files
    embedding_path = 'example_data/';
    
    %path to *Frames.txt files containing  
    led_file_path = 'example_data/';
    
    %path to training set file
    training_set_path = 'example_data/DI_training_set.mat';
    
    %     %path to example movies
    %     example_movie_path = '/Volumes/chrimson/video_library/';
    %
    %     %path to example movie information .mat file
    %     example_movie_info_path = '/Volumes/chrimson/video_library/exampleMovieInformation.mat';
        
    %region boundaries for the behavioral space map
    xRange = [-105 105];
    
    %number of points in each dimension of the map
    numPoints = 101;
    
    %image sampling rate (in Hz)
    Fs = 100;
    
    %use entropy window or not 
    use_entropy_window = false;
    
    %Time width for entropy estimation (in frames)
    entropyTimeWidth = 20;
    
    %Number of bins for entropy estimation
    entropyNumBins = 51;
    
    %Starting time for entropy plot
    entropyStartTime = -22.5;
    
    %Threshold for density calculation (density < 1e-7 = 0)
    densityThreshold = 1e-7;
    
    %statistical signficance level (after correcting for multiple
    %hypotheses)
    sigAlpha = .05;
    
    %disk size for morphological closing of identified regions
    closeRegionSize = 2;
    
    %disk size for morphological dilation of identified regions
    dilateRegionSize = 1;
    
    %time spacing of points in stimulus-triggered plot
    triggered_plot_time_spacing = 20;
    
    %time width for averaging in stimulus-triggered plot
    triggered_plot_time_width = 50;
    
    %Minimum size for a region (in pixels)
    minRegionSize = 10;
    
    %statistical significance for time series plot error bars
    numSigma = 1;
    
    %number of time bins for time series plots
    timeBins = 300;
    
    %number of bins for calculating time series values
    yBins = 100;
    
    %window for making on/off comparisons
    off_window = [30 45];
    
    %smoothing length scale for kernel density estimation
    sigma = 1.5;
    
    %maximum number of LED cycles to use (negative number implies use all)
    maxNumCycles = -1;    
    
    %use 'session' to analyze a single filming session or 'strain' to 
    %analyze all sessions of the strain (currently unused in this code)
    session_or_strain = 'session';
    
    %Number of bootstrap simulations for time series plots
    bootstrapNum = 1000;
    
    %Smoothing size for time series plots (in seconds)
    timeBinSize = 3;
    
    %Fractional confidence bound for time series plots
    confBound = .67;
    
    %display plots?
    plotsOn = true;
    
    %ensure that the entropy curve is calculated
    force_calculate_entropy = false;
    
    %set all data sets to be experimentals
    no_controls = false;
    
    
    %%%%%%%%%%%%%%%%%%
    
    
    if ~isfield(parameters,'no_controls') || isempty(parameters.no_controls)
        parameters.no_controls = no_controls;
    end
    
    if ~isfield(parameters,'force_calculate_entropy') || isempty(parameters.force_calculate_entropy)
        parameters.force_calculate_entropy = force_calculate_entropy;
    end
    
    if ~isfield(parameters,'use_entropy_window') || isempty(parameters.use_entropy_window)
        parameters.use_entropy_window = use_entropy_window;
    end
    
    if ~isfield(parameters,'session_or_strain') || isempty(parameters.session_or_strain)
        parameters.session_or_strain = session_or_strain;
    end
    
    if ~isfield(parameters,'embedding_path') || isempty(parameters.embedding_path)
        parameters.embedding_path = embedding_path;
    end
    
    if ~isfield(parameters,'led_file_path') || isempty(parameters.led_file_path)
        parameters.led_file_path = led_file_path;
    end
    
    if ~isfield(parameters,'training_set_path') || isempty(parameters.training_set_path)
        parameters.training_set_path = training_set_path;
    end
    
    if ~isfield(parameters,'xRange') || isempty(parameters.xRange)
        parameters.xRange = xRange;
    end
    
    if ~isfield(parameters,'numPoints') || isempty(parameters.numPoints)
        parameters.numPoints = numPoints;
    end
    
    if ~isfield(parameters,'Fs') || isempty(parameters.Fs)
        parameters.Fs = Fs;
    end
    
    if ~isfield(parameters,'entropyTimeWidth') || isempty(parameters.entropyTimeWidth)
        parameters.entropyTimeWidth = entropyTimeWidth;
    end
    
    if ~isfield(parameters,'entropyNumBins') || isempty(parameters.entropyNumBins)
        parameters.entropyNumBins = entropyNumBins;
    end
    
    if ~isfield(parameters,'entropyStartTime') || isempty(parameters.entropyStartTime)
        parameters.entropyStartTime = entropyStartTime;
    end
    
    if ~isfield(parameters,'off_window') || isempty(parameters.off_window)
        parameters.off_window = off_window;
    end
    
    if ~isfield(parameters,'densityThreshold') || isempty(parameters.densityThreshold)
        parameters.densityThreshold = densityThreshold;
    end
    
    if ~isfield(parameters,'sigAlpha') || isempty(parameters.sigAlpha)
        parameters.sigAlpha = sigAlpha;
    end
    
    if ~isfield(parameters,'closeRegionSize') || isempty(parameters.closeRegionSize)
        parameters.closeRegionSize = closeRegionSize;
    end
    
    if ~isfield(parameters,'triggered_plot_time_spacing') || isempty(parameters.triggered_plot_time_spacing)
        parameters.triggered_plot_time_spacing = triggered_plot_time_spacing;
    end
    
    if ~isfield(parameters,'triggered_plot_time_width') || isempty(parameters.triggered_plot_time_width)
        parameters.triggered_plot_time_width = triggered_plot_time_width;
    end
    
    if ~isfield(parameters,'minRegionSize') || isempty(parameters.minRegionSize)
        parameters.minRegionSize = minRegionSize;
    end
    
    if ~isfield(parameters,'numSigma') || isempty(parameters.numSigma)
        parameters.numSigma = numSigma;
    end
    
    if ~isfield(parameters,'timeBins') || isempty(parameters.timeBins)
        parameters.timeBins = timeBins;
    end
    
    if ~isfield(parameters,'yBins') || isempty(parameters.yBins)
        parameters.yBins = yBins;
    end
    
    if ~isfield(parameters,'sigma') || isempty(parameters.sigma)
        parameters.sigma = sigma;
    end    
    
    if ~isfield(parameters,'dilateRegionSize') || isempty(parameters.dilateRegionSize)
        parameters.dilateRegionSize = dilateRegionSize;
    end
    
    %     if ~isfield(parameters,'example_movie_path') || isempty(parameters.example_movie_path)
    %         parameters.example_movie_path = example_movie_path;
    %     end
    %
    %     if ~isfield(parameters,'example_movie_info_path') || isempty(parameters.example_movie_info_path)
    %         parameters.example_movie_info_path = example_movie_info_path;
    %     end
    
    if ~isfield(parameters,'maxNumCycles') || isempty(parameters.maxNumCycles)
        parameters.maxNumCycles = maxNumCycles;
    end
    
    if ~isfield(parameters,'bootstrapNum') || isempty(parameters.bootstrapNum)
        parameters.bootstrapNum = bootstrapNum;
    end
    
    if ~isfield(parameters,'timeBinSize') || isempty(parameters.timeBinSize)
        parameters.timeBinSize = timeBinSize;
    end
    
    if ~isfield(parameters,'confBound') || isempty(parameters.confBound)
        parameters.confBound = confBound;
    end
    
    if ~isfield(parameters,'plotsOn') || isempty(parameters.plotsOn)
        parameters.plotsOn = plotsOn;
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    