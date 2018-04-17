function xVals = labelEntropyPeak(ts,ent_exp,ent_control)


    if ~isempty(ent_control)
        plot(ts,ent_exp,'r-',ts,ent_control,'b-')
    else
        plot(ts,ent_exp,'r-')
    end
    
    hold on
    xlim([-10 25]);
    w = ylim;
    legend({'Exp','Control'},'fontsize',14,'fontweight','bold')
    xlabel('Time (s)','fontsize',14,'fontweight','bold')
    ylabel('Entropy (bits)','fontsize',14,'fontweight','bold')
    set(gca,'fontsize',12,'fontweight','bold')
    
    test = true;
    
    while test
        
        title('Label Start Position','fontsize',16,'fontweight','bold');
        q = -1;
        while q ~= 1
            [x1,~,q] = ginput(1);
        end
        plot(x1.*[1 1],[-10 25],'k--');
        ylim(w)
        
        title('Label End Position','fontsize',16,'fontweight','bold');
        q = -1;
        while q ~= 1
            [x2,~,q] = ginput(1);
        end
        plot(x2.*[1 1],[-10 25],'k--');
        ylim(w)
        
        x = [x1 x2 x2 x1 x1];
        y = w([1 1 2 2 1]);
        fill(x,y,[0 1 0],'facealpha',.2,'edgealpha',0);
        
        title('OK? (y/n)','fontsize',16,'fontweight','bold')
        q = -1;
        while (q ~= 121) && (q ~= 110)
            [~,~,q] = ginput(1);
            if isempty(q)
                q = 9999;
            end
        end
        
        if q == 121
            test = false;
        end
        
        if test
            clf
            if ~isempty(ent_control)
                plot(ts,ent_exp,'r-',ts,ent_control,'b-')
            else
                plot(ts,ent_exp,'r-')
            end
            
            hold on
            xlim([-10 25]);
            w = ylim;
            legend({'Exp','Control'},'fontsize',14,'fontweight','bold')
            xlabel('Time (s)','fontsize',14,'fontweight','bold')
            ylabel('Entropy (bits)','fontsize',14,'fontweight','bold')
            set(gca,'fontsize',12,'fontweight','bold')
        end
        
    end
    drawnow
    
    xVals = [x1 x2];