clear all;close all;

addpath(path,'../../matlab_tools/m_map/');
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

method = 'mean';
mat_dat = 'C:\Users\cheny\Desktop\EcoHAB\Cruise_data\FWC_kb_dat_0006.mat';

for check_time = datenum(2005,5,1):1:datenum(2023,12,1)

    t_min = check_time;
    t_max = check_time;

    try
        [kb_lon,kb_lat, kb_dep, kb_dat, kb_dep_all, kb_dat_all] = get_FWC_kb(t_min,t_max,mat_dat,method);
    catch
        continue;
    end

    kb_dat(kb_dat<10) = NaN;

    if(sum(~isnan(kb_dat)))
        figure(1);
        m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
        m_scatter(kb_lon,kb_lat,50,log10(kb_dat),'filled');
        caxis([3 7]);
        hold on;
        m_gshhs_h('patch',[.6 .6 .6]);
        hold on;
        m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
        colorbar;
        set(gcf,'color', [1 1 1]);
        xlabel('Longitude','fontsize',9);
        ylabel('Latitude','fontsize',9);
        title(datestr(check_time,'yyyy-mmm-dd'));
        hold off;
        pause
    end
end