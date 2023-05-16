function [ax] = plot_var_flow2(fig_path,t1,v1,t2,v2,var_name,name,psize,xlab,ylab,range)
%   Detailed explanation goes here

    left = 0.07;
    bottom = 0.1;
    width = 1/psize(2)*0.75;
    height = 1/psize(1)*0.7;
    blk_h = 0.055;
    blk_v = 0.055;
    
    k=0;
    k2=0;
    for i=1:size(v2,2)
        if(mod(k,psize(1)*psize(2))==0)
            close(gcf);
            figure('Units','pixels','position',[10 10 900 1200])
            set(gcf,'color',[1 1 1]);
            k=0;
        end
        k=k+1;
        k2=k2+1;
        row = ceil(k/psize(2));
        col = mod(k,psize(2));
        if(col==0)
            col=psize(2);
        end
        ax(k) = axes('units','normalized','position',[left+(col-1)*(width+blk_h) bottom+(psize(1)-row)*(height+blk_v) width height]);   
        x = v1{i};
        y = v2{i};
        scatter(x,y,'k','filled');
        hold on;
        x = c1{i};
        y = c2{i};
        plot(x,y,'linewidth',1.5);
        legend({'Data','AMLE','MLE','LAD'},'location','southeast');
        try
            text(0.83,0.88,var_name{i},'units','normalized','fontsize',11,'fontweight','bold');
        catch
            text(0.83,0.88,'No Data','units','normalized','fontsize',11,'fontweight','bold');
        end
        
%         if(col>1&&row<psize(1))
%             set(gca,'xticklabel',[]);
%             set(gca,'yticklabel',[]);
%         elseif(col==1&&row<psize(1))
%             set(gca,'xticklabel',[]);
%             ylabel(ylab);
%         elseif(col>1&&row==psize(1))
%             set(gca,'yticklabel',[]);
%             xlabel(xlab);
%         else
%             ylabel(ylab);
%             xlabel(xlab);
%         end
        
        ylabel(ylab);
        xlabel(xlab);
        
        if(~isempty(range))
            axis(range);
        end
        set(gca,'fontsize',10);
        if(mod(k,psize(1)*psize(2))==0||i==size(v2,2))   
            outname = [fig_path name '_' num2str(ceil(k2/(psize(1)*psize(2))))];
            print(gcf,outname,'-dpng'); 
        end
    end
end