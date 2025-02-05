clear all;close all;

addpath(path,'../../matlab_tools/m_map/');
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

method = 'mean';
k=0;
figure('Units','pixels','Position',[50 50 920 440]);

tiledlayout(1,2,"TileSpacing","compact","Padding","compact");
for check_time = [datenum(2022,10,07) datenum(2022,10,17)]
    k=k+1;
    t_min = check_time;
    t_max = check_time;

    try
        [kb_lon,kb_lat, kb_dep, kb_dat, kb_dep_all, kb_dat_all] = get_FWC_kb(t_min,t_max,method);
    catch
        continue;
    end

    kb_dat(kb_dat<10) = NaN;

    nexttile;
    m_proj('Mercator','lat',[26.5 28.0],'long',[-83.6 -81.9]);
    m_scatter(kb_lon,kb_lat,200,kb_dat,'^','filled');
    hold on;
    m_gshhs_f('patch',[.8 .8 .8]);
    hold on
    m_grid('linestyle','none','linewidth',1,'fontsize',14);
    t = colorbar;
    t.Units = 'normalized';
    t.Location = 'north';
    set(t,'FontSize',14);
    if(k==1)
        t.Position = [0.12 0.25 0.15 0.04];
        ylabel('Latitude','fontsize',14);
    elseif(k==2)
        t.Position = [0.62 0.25 0.15 0.04];
    end
    xlabel(t,'cells/L')
    set(gcf,'color', [1 1 1]);
    
    text('Units','normalized','Position',[0.78 0.95],'String',datestr(check_time,'mmm-dd'),'FontSize',18);

    hold off;
end

exportgraphics(gcf,'fwc_data.png','Resolution',300);