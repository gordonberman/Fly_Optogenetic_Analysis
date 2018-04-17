function makeTimeSeriesPlot(ts,vals,LEDs,isControl,...
                    timeBinSize,bootstrapNum,confBound) 

                        
    if nargin < 5 || isempty(timeBinSize)
        timeBinSize = 5;
    end
                        
    if nargin < 6 || isempty(bootstrapNum)
        bootstrapNum = 1000;
    end
  
    if nargin < 7 || isempty(confBound) 
        confBound = .67;
    end
    

    confBound = 1 - confBound;
    dt = 1/100;
    cm = [200 200 200]./255;
    alphaValue = .4;
    
    ts = ts .* dt;
    N = length(LEDs);
    dt2 = ts(2) - ts(1);
    alpha = 1 ./ (timeBinSize/dt2);
    
    L = length(ts);
    timeSeries = zeros(N,L);
    for i=1:N       
        x = vals(i,:);
        x(isnan(x) | isinf(x)) = 0;
        timeSeries(i,:) = causalfilterdata(x,alpha);
    end
    
    
    
    if sum(isControl) > 0
        meanControl = mean(timeSeries(isControl,:));
        a = timeSeries(isControl,:);
        while length(a(:,1)) < 3
            b = [a;a];
            a = b;
        end
        controlBounds = bootci(bootstrapNum,{@mean,a},'alpha',confBound);
    end
        
    meanExp = mean(timeSeries(~isControl,:));
    a = timeSeries(~isControl,:);
    while length(a(:,1)) < 3
        b = [a;a];
        a = b;
    end
    expBounds = bootci(bootstrapNum,{@mean,a},'alpha',confBound);
    
    if sum(isControl) > 0
        maxVal = max(max(controlBounds(:)),max(expBounds(:)));
    else
        maxVal = max(expBounds(:));
    end
    minVal = 0;
    
    
    figure
    CC = bwconncomp(LEDs{1});
    hold on
    for i=1:CC.NumObjects
        x = CC.PixelIdxList{i}(1)*dt./60;
        y = CC.PixelIdxList{i}(end)*dt./60;
        rectangle('Position',[x minVal y-x+dt maxVal-minVal],'facecolor',cm,'edgecolor',cm);
    end
    
    x = ts./60;
    x(1) = 0;
    X = [x fliplr(x) x(1)];
    Yexp = [expBounds(1,:) fliplr(expBounds(2,:)) expBounds(1,1)];
    if sum(isControl) > 0
        Ycontrol = [controlBounds(1,:) fliplr(controlBounds(2,:)) controlBounds(1,1)];
        fill(X,Ycontrol,'b','facealpha',alphaValue,'edgealpha',0);
    end   
    
    fill(X,Yexp,'r','facealpha',alphaValue,'edgealpha',0)
    
    if sum(isControl) > 0
        plot(x,meanControl,'b-','linewidth',2);
    end
    plot(x,meanExp,'r-','linewidth',2);
    
    xlim([0 max(ts)./60]);
    ylim([0 maxVal]);
    
    xlabel('Time (min)','fontsize',16,'fontweight','bold')
    ylabel('Region Density','fontsize',16,'fontweight','bold')
    legend({'Controls','Experimentals'},'fontsize',16,'fontweight','bold');
    set(gca,'fontsize',14,'fontweight','bold')
