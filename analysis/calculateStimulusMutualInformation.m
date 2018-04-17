function MI_data = calculateStimulusMutualInformation(densities,assignments,averageDensity,numRegions,xx)


    if iscell(assignments)
        assignments = cell2mat(assignments);
    end
    
    
    if iscell(densities)
        L = length(densities);
        densities = cell2mat(reshape(densities,[1 1 L]));
    end
    

    numPoints = length(xx);
    dx = xx(2) - xx(1);

    numRegionAssignments = zeros(numRegions+1,1);
    averageRegions = zeros(numPoints,numPoints,numRegions+1);
    for j=1:numRegions+1
        idx = assignments==j-1;
        averageRegions(:,:,j) = sum(densities(:,:,idx),3);
        numRegionAssignments(j) = sum(idx);
    end
    
    for i=1:numRegions+1
        averageRegions(:,:,i) = averageRegions(:,:,i) / numRegionAssignments(i);
    end
    region_probs = numRegionAssignments / sum(numRegionAssignments);
    
    
    MI_data.MI = 0;
    for i=1:numRegions+1
        z = averageRegions(:,:,i).*log2(averageRegions(:,:,i)./averageDensity);
        z(isnan(z) | isinf(z)) = 0;
        MI_data.MI = MI_data.MI + sum(z(:))*dx^2*region_probs(i);
    end
    
    MI_data.region_probs = region_probs;
    MI_data.averageRegions = averageRegions;
    