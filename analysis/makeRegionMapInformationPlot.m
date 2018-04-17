function makeRegionMapInformationPlot(MI_data)

    load('saved_colormaps.mat','cc');
    
    N = MI_data.numRegions;
    M = ceil((N+3)/2);
    
    cs = 'brgmck';
    
    subplot(M,2,1)   
    regionLabels = cell(N,1);
    hold on
    w = zeros(N,1);
    plot(MI_data.region_boundary{1}(:,2),MI_data.region_boundary{1}(:,1),'k-','linewidth',2)
    for i=1:N
        B = bwboundaries(MI_data.originalRegionMap == i);
        w(i)= fill(B{1}(:,2),B{1}(:,1),cs(mod(i-1,6)+1),'facealpha',.9,'edgecolor',cs(mod(i-1,6)+1),'edgealpha',1);
        regionLabels{i} = ['Region #' num2str(i)];
    end
    legend(w,regionLabels,'fontsize',12,'fontweight','bold','location','eastoutside')
    strain_name = MI_data.sessionName;
    idx = find(strain_name == '_',1,'first');
    if ~isempty(idx)
        strain_name = strain_name(1:idx-1);
    end
    title(strain_name,'fontsize',16,'fontweight','bold')
    
    axis equal tight xy off
    
    
    for i=1:(N+1)
        subplot(M,2,i+2)
        imagesc(MI_data.average_prior_regions_exp(:,:,i))
        axis equal tight off xy
        hold on
        colormap(cc)
        plot(MI_data.region_boundary{1}(:,2),MI_data.region_boundary{1}(:,1),'k-','linewidth',2)
        caxis([0 5e-4])
        colorbar
        set(gca,'fontsize',12,'fontweight','bold')
        q = num2str(round(100*MI_data.prob_after_region_assignments_exp(i))/100);
        title(['p( x | R_' num2str(i-1) ' ),  p( R_' num2str(i-1) ' ) = ' q ],'fontsize',14,'fontweight','bold')   
    end
    
    
    partial_infos = zeros(N+1,1);
    xx = MI_data.xx;
    for i=1:(N+1)
        z = MI_data.average_prior_regions_exp(:,:,i).*log2(MI_data.average_prior_regions_exp(:,:,i)./MI_data.density_exp);
        z(isinf(z) | isnan(z)) = 0;
        partial_infos(i) = sum(z(:))*(xx(2)-xx(1))^2*MI_data.prob_after_region_assignments_exp(i);
    end
    
    if N > 1
        partial_infos_no_zeros = zeros(N,1);
        xx = MI_data.xx;
        p = MI_data.prob_after_region_assignments_exp_no_zero;
        p(1) = 0;
        p = p ./ sum(p);
        for i=2:(N+1)
            z = MI_data.average_prior_regions_exp(:,:,i).*log2(MI_data.average_prior_regions_exp(:,:,i)./MI_data.density_exp);
            z(isinf(z) | isnan(z)) = 0;
            partial_infos_no_zeros(i-1) = sum(z(:))*(xx(2)-xx(1))^2*p(i);
        end
    end
    
    
    subplot(M,2,2)
    if N > 1
        
        plot(0:N,partial_infos,'bo','markerfacecolor','b');
        hold on
        plot(1:N,partial_infos_no_zeros,'rs','markerfacecolor','r');
        set(gca,'fontsize',12,'fontweight','bold','xtick',0:N)
        xlim([-.5 N+.5])
        xlabel('Region #','fontsize',14,'fontweight','bold')
        ylabel('Partial Mutual Information (bits)','fontsize',14,'fontweight','bold')
        legend({'Including Insignificant Regions','Without Insignificant Regions'},...
            'fontsize',12,'fontweight','bold','location','southoutside')
        q = ylim;
        ylim([0 q(2)]);
        
        q = num2str(round(MI_data.MI_estimate_exp.MI_estimate*100)/100);
        
        q2 = num2str(round(100*MI_data.MI_estimate_exp_no_zero.MI_estimate)/100);
        
        title({['Mutual Information = ' q ' bits'],['(' q2 ' bits Between Regions)']},'fontsize',16,'fontweight','bold') 
        
                
    else
        
        plot(0:N,partial_infos,'bo','markerfacecolor','b');
        set(gca,'fontsize',12,'fontweight','bold','xtick',0:N)
        xlim([-.5 N+.5])
        xlabel('Region #','fontsize',14,'fontweight','bold')
        ylabel('Partial Mutual Information (bits)','fontsize',14,'fontweight','bold')
        q = ylim;
        ylim([0 q(2)]);
        
        q = num2str(round(MI_data.MI_estimate_exp.MI_estimate*100)/100);
        q = q(1:min(4,length(q)));
        
        title(['Mutual Information = ' q ' bits'],'fontsize',16,'fontweight','bold') 
                        
    end
    
    
    
    