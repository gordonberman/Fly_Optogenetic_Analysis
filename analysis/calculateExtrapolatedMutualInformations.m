function outputData = calculateExtrapolatedMutualInformations(assignments,densities,numRegions,maxNumPartitions,partitionReplicates,xx)


    ms = (2:maxNumPartitions)';
    outputData.ms = ms;
    L = length(assignments);
    
    MI_values = cell(length(ms),1);
    MI_means = zeros(length(ms),1);
    MI_vars = zeros(length(ms),1);
    for i=1:length(ms)
        q = zeros(partitionReplicates,ms(i));
        m = ms(i);
        parfor j=1:partitionReplicates
            [~,idx] = sort(randperm(L));
            z = assignments;
            d = densities;
            M = ceil(L/m);
            rr = zeros(1,m);
            for k=1:m
                w = (1:M) + (k-1)*M;
                w = idx(w(w<=L));
                x = calculateStimulusMutualInformation(d(:,:,w),...
                    z(w),mean(d(:,:,w),3),numRegions,xx);
                rr(k) = x.MI;
            end
            q(j,:) = rr;
        end
        MI_values{i} = q;
        MI_means(i) = mean(q(:));
        MI_vars(i) = var(q(:));
    end
    
    
    
    p_mean = polyfit(L./(ms),MI_means,1);
    p_var = polyfit(log(ms),log(MI_vars),1);
    
    outputData.MI_values = MI_values;
    outputData.MI_means = MI_means;
    outputData.MI_vars = MI_vars;
    outputData.p_mean = p_mean;
    outputData.p_var = p_var;
    
    outputData.MI_estimate = p_mean(2);
    outputData.MI_std = sqrt(exp(p_var(2)));
    outputData.L = L;