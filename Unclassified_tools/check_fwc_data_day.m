clear all;close all;

addpath(path,'../../../matlab_tools/m_map/');
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

method = 'mean';

for plot_date = datenum(2022,9,10):datenum(2022,10,31)

    t_min = plot_date;
    t_max = plot_date+1;

    try
        [kb_lon,kb_lat, kb_dep, kb_dat, kb_dep_all, kb_dat_all] = get_FWC_kb(t_min,t_max,method);
    catch
        continue;
    end  

    kb_dat(kb_dat<10) = NaN;

    if(sum(~isnan(kb_dat))==0)
        continue;
    end

    figure(1);
    m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
    m_scatter(kb_lon,kb_lat,50,log10(kb_dat),'filled');
    hold on;
    m_gshhs_h('patch',[.6 .6 .6]);
    hold on;
    m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
    colorbar;
    set(gcf,'color', [1 1 1]);
    xlabel('Longitude','fontsize',9);
    ylabel('Latitude','fontsize',9);
    title(datestr(plot_date));
    

    hold off;
    pause
end