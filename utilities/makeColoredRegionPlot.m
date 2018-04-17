function makeColoredRegionPlot(images,xx,alphaValue,cMap,dilateSize)


    if nargin < 3 || isempty(alphaValue)
        alphaValue = .4;
    end

    if nargin < 4 || isempty(cMap)
        cMap = jet;
    end
    
    if nargin < 5 || isempty(dilateSize)
        dilateSize = 2;
    end

    N = length(images(1,1,:));
    
    cx = fit((1:64)',cMap(:,1),'linearinterp');
    cy = fit((1:64)',cMap(:,2),'linearinterp');
    cz = fit((1:64)',cMap(:,3),'linearinterp');
    
    hold on;
    for i=1:N
        
        z = 63*((i-1)/(N-1)) + 1;
        
        if dilateSize > 0
           A = imclose(images(:,:,i)>0,strel('disk',dilateSize));
        else
           A = images(:,:,i) > 0;
        end
        
        c = rectify([cx(z) cy(z) cz(z)]);
        c(c > 64) = 64;
        
        B = bwboundaries(A);
        for j=1:length(B)
            fill(xx(B{j}(:,2)),xx(B{j}(:,1)),c,...
                'facealpha',alphaValue,'edgealpha',0);
        end
        
    end
    
    axis equal tight off xy;