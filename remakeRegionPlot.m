function stimulusPlots = remakeRegionPlot(outputStats,regionNum,omitTimeSeriesPlots)
%remakes output plots from any number of regions
%Outputs:
%   stimulusPlots -> cell array of outputs from makeStimulusTriggeredPlots
%
%Inputs:
%   outputStats -> struct containing output from analyzeSessions
%   regionNum -> array of regions for which to remake plots


    addpath('utilities');
    addpath('subroutines');
    load('saved_colormaps');
    
    regions = outputStats.regions;
    numSigma = outputStats.parameters.numSigma;
    timeBins = outputStats.parameters.timeBins;
    timeBinSize = outputStats.parameters.timeBinSize;
    bootstrapNum = outputStats.parameters.bootstrapNum;
    confBound = outputStats.parameters.confBound;
    yBins = outputStats.parameters.yBins;
    %BB = bwboundaries(outputStats.density > outputStats.densityThreshold);
    load('DI_training_set.mat','B');
    BB = B;
    vals = outputStats.vals;
    ts = outputStats.ts;
    isControl = outputStats.isControl;
    startTime = outputStats.startTime;
    xx = outputStats.xx;
    dx = xx(2) - xx(1);
    endTime = outputStats.endTime;
    LEDs = outputStats.LEDs;
    numRegions = outputStats.numRegions;
    windows = outputStats.windows;
    window_densities = outputStats.window_densities;
    LED_spacing = outputStats.LED_spacing;
    
    if nargin < 2 || isempty(regionNum)
        regionNum = 1:length(regions);
    end
    
    if nargin < 3 || isempty(omitTimeSeriesPlots)
        omitTimeSeriesPlots = false;
    end
    
    stimulusPlots = cell(sum(numRegions),1);
    cSums = [0; cumsum(numRegions)];
    count = 1;
    
    for k=1:length(numRegions)
        
        for i=1:numRegions(k)
            
            figure
            r = cSums(k) + i;
            
            subplot(1,3,2)
            currentVals = vals(:,:,r);
            stimulusPlots{count} = makeStimulusTriggeredPlots(ts,currentVals,...
                LEDs,isControl,numSigma,timeBins,yBins,startTime,endTime,LED_spacing);
            count = count + 1;
            
            subplot(1,3,3)
            M = outputStats.regionImages(:,:,r);
            M = sqrt(M.^2 ./ sum(M(:).^2));
            imagesc(M);
            colormap(cc)
            
            caxis([.01 .03])
            axis equal tight off xy
            freezeColors;
            
            subplot(1,3,1)
            
            imagesc(xx,xx,regions{k} == i);
            
            hold on
            axis equal tight off xy
            
            B = bwboundaries(regions{k} == i);
            for j=1:length(B)
                plot(xx(B{j}(:,2)),xx(B{j}(:,1)),'k-','linewidth',2);
            end
            
            for j=1:length(BB)
                plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',3)
            end
            
            if k == 1
                title(['+ Region #' num2str(i)],'fontsize',16,'fontweight','bold')
                colormap(cc)
            else
                title(['- Region #' num2str(i)],'fontsize',16,'fontweight','bold')
                colormap(cc3)
            end
            
            freezeColors;
            
            
            drawnow
            
            if omitTimeSeriesPlots
                close(gcf)
            end
            
            if ~omitTimeSeriesPlots
                
                makeTimeSeriesPlot(ts,currentVals,LEDs,isControl,timeBinSize,bootstrapNum,confBound);
                if k == 1
                    title(['+ Region #' num2str(i)],'fontsize',16,'fontweight','bold')
                    colormap(cc)
                else
                    title(['- Region #' num2str(i)],'fontsize',16,'fontweight','bold')
                    colormap(cc3)
                end
                
                drawnow
                
            end
            
        end
        
    end
    
    
    figure
    a = round(windows.*100)/100;
    
    
    
    for i=1:4
        if sum(isControl) > 0
            subplot(2,4,i)
        else
            subplot(1,4,i)
        end
        imagesc(xx,xx,window_densities(:,:,i))
        axis equal tight off xy
        hold on;
        for j=1:length(BB)
            plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',2)
        end
        title([num2str(a(i,1)) 's to ' num2str(a(i,2)) 's'],'fontsize',14,'fontweight','bold')
        colormap(cc)
        caxis([0 4e-4])
    end
    
    if sum(isControl) > 0
        for k=1:4
            i = k+4;
            subplot(2,4,i)
            imagesc(xx,xx,window_densities(:,:,i))
            axis equal tight off xy
            hold on;
            for j=1:length(BB)
                plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',2)
            end
            title([num2str(a(k,1)) 's to ' num2str(a(k,2)) 's (control)'],'fontsize',14,'fontweight','bold')
            colormap(cc)
            caxis([0 4e-4])
        end
    end
    
    if omitTimeSeriesPlots
        close(gcf)
    end
    
    
    colorVals = 'brkmgc';
    segs = {'-','--','-.'};
    
    if sum(isControl) > 0 && ~isempty(numRegions)
        if numRegions(1) > 0
            
            if sum(isControl) > 0
                figure
                hold on
                hs = zeros(numRegions(1),1);
                for i=1:numRegions(1)
                    tVals = stimulusPlots{i}.ts;
                    m = stimulusPlots{i}.m;
                    s = stimulusPlots{i}.s;
                    m2 = stimulusPlots{i}.m2;
                    s2 = stimulusPlots{i}.s2;
                    
                    m = m - m2;
                    s = sqrt(s.^2 + s2.^2);
                    
                    X = [tVals(1:end-1) fliplr(tVals(1:end-1)) tVals(1)];
                    Y = [m(1:end-1)-s(1:end-1); flipud(m(1:end-1)+s(1:end-1));m(1)-s(1)];
                    
                    hs(i) = plot(tVals,m,[colorVals(mod(i-1,6)+1) segs{floor((i-1)/5)+1}],'linewidth',2);
                    fill(X,Y,colorVals(mod(i-1,6)+1),'edgealpha',0,'facealpha',.5)
                    
                end
                yvals = ylim;
                plot([0 0],yvals,'k--','linewidth',1)
                plot([endTime endTime],yvals,'k--','linewidth',1)
                plot(stimulusPlots{i}.ts,zeros(size(stimulusPlots{i}.ts)),'k--','linewidth',1)
                ylim(yvals);
                set(gca,'fontsize',14,'fontweight','bold')
                xlabel('Time From Stimulus (s)','fontsize',16,'fontweight','bold')
                ylabel('Change in Region Density Above Controls','fontsize',16,'fontweight','bold')
                a = repmat({'Region #'},numRegions(1),1);
                for i=1:length(a)
                    a{i} = [a{i} num2str(i)];
                end
                legend(hs,a,'fontsize',16,'fontweight','bold','location','northwest')
            end
            
            figure
            hold on
            hs = zeros(numRegions(1),1);
            for i=1:numRegions(1)
                tVals = stimulusPlots{i}.ts;
                m = stimulusPlots{i}.m;
                s = stimulusPlots{i}.s;
                
                X = [tVals(1:end-1) fliplr(tVals(1:end-1)) tVals(1)];
                Y = [m(1:end-1)-s(1:end-1); flipud(m(1:end-1)+s(1:end-1));m(1)-s(1)];
                
                hs(i) = plot(tVals,m,[colorVals(mod(i-1,6)+1) segs{floor((i-1)/5)+1}],'linewidth',2);
                fill(X,Y,colorVals(mod(i-1,6)+1),'edgealpha',0,'facealpha',.5)
                
            end
            yvals = ylim;
            plot([0 0],yvals,'k--','linewidth',1)
            plot([endTime endTime],yvals,'k--','linewidth',1)
            plot(stimulusPlots{i}.ts,zeros(size(stimulusPlots{i}.ts)),'k--','linewidth',1)
            ylim(yvals);
            set(gca,'fontsize',14,'fontweight','bold')
            xlabel('Time From Stimulus (s)','fontsize',16,'fontweight','bold')
            ylabel('Total Region Density','fontsize',16,'fontweight','bold')
            a = repmat({'Region #'},numRegions(1),1);
            for i=1:length(a)
                a{i} = [a{i} num2str(i)];
            end
            legend(hs,a,'fontsize',16,'fontweight','bold','location','northwest')
            
            
            figure
            hold on
            hs = zeros(numRegions(1),1);
            for i=1:numRegions(1)
                
                regionSize = sum(sum(regions{1}==i))*dx^2;
                
                tVals = stimulusPlots{i}.ts;
                m = stimulusPlots{i}.m ./ regionSize;
                s = stimulusPlots{i}.s ./ regionSize;
                        
                X = [tVals(1:end-1) fliplr(tVals(1:end-1)) tVals(1)];
                Y = [m(1:end-1)-s(1:end-1); flipud(m(1:end-1)+s(1:end-1));m(1)-s(1)];
                
                hs(i) = plot(tVals,m,[colorVals(mod(i-1,6)+1) segs{floor((i-1)/5)+1}],'linewidth',2);
                fill(X,Y,colorVals(mod(i-1,6)+1),'edgealpha',0,'facealpha',.5)
                
            end
            yvals = ylim;
            plot([0 0],yvals,'k--','linewidth',1)
            plot([endTime endTime],yvals,'k--','linewidth',1)
            plot(stimulusPlots{i}.ts,zeros(size(stimulusPlots{i}.ts)),'k--','linewidth',1)
            ylim(yvals);
            set(gca,'fontsize',14,'fontweight','bold')
            xlabel('Time From Stimulus (s)','fontsize',16,'fontweight','bold')
            ylabel('Average Region Density','fontsize',16,'fontweight','bold')
            a = repmat({'Region #'},numRegions(1),1);
            for i=1:length(a)
                a{i} = [a{i} num2str(i)];
            end
            legend(hs,a,'fontsize',16,'fontweight','bold','location','northwest')
            
            
            figure
            
            toPlot = setdiff(unique(regions{1}),0);
            L = length(toPlot);
            N = ceil(sqrt(L));
            M = ceil(L/N);
            for i=1:L
                subplot(M,N,i)
                
                Q = zeros(size(regions{1}));
                Q(regions{1} == i) = i;
                
                imagesc(xx,xx,Q);
                
                hold on
                axis equal tight off xy
                
                B = bwboundaries(regions{1} == i);
                for j=1:length(B)
                    plot(xx(B{j}(:,2)),xx(B{j}(:,1)),'k-','linewidth',2);
                end
                
                for j=1:length(BB)
                    plot(xx(BB{j}(:,2)),xx(BB{j}(:,1)),'k-','linewidth',3)
                end
                
                title(['Region #' num2str(i)],'fontsize',16,'fontweight','bold')
                colormap(cc6)
                caxis([0 8])
                
            end
        end
    end
    