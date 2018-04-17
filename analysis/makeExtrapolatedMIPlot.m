function makeExtrapolatedMIPlot(MI_output)

    L = MI_output.L;
    ms = MI_output.ms;
    MI_means = MI_output.MI_means;
    MI_vars = MI_output.MI_vars;
    p_mean = MI_output.p_mean;
    p_var = MI_output.p_var;
    
    outputData.MI_estimate = p_mean(2);
    outputData.MI_std = sqrt(exp(p_var(2)));

    figure
    subplot(1,2,1)
    plot(L./(ms),MI_means,'ks','markerfacecolor','k')
    hold on
    plot(L./(ms),polyval(p_mean,L./(ms)),'r-','linewidth',1)
    errorbar(0,p_mean(2),outputData.MI_std,'ro','markerfacecolor','r')
    xlabel('1/N')
    ylabel('Mutual Information (bits)')
    set(gca,'fontsize',12,'fontweight','bold')
    
    
    subplot(1,2,2)
    plot(log(ms),log(MI_vars),'ks','markerfacecolor','k')
    hold on
    plot(log(ms),polyval(p_var,log(ms)),'r-','linewidth',1)
    plot(0,p_var(2),'ro','markerfacecolor','r')
    xlabel('Log(m)')
    ylabel('Log(Var(Mutual Information) (bits^2))')
    set(gca,'fontsize',12,'fontweight','bold')