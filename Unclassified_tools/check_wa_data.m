clear all;close all;

addpath(path,'../../matlab_tools/m_map/');
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

method = 'mean';

load("WA_chla_dat.mat");

for check_time = datenum(2022,10,1):datenum(2022,12,30)

    pos = find(s_date>=check_time&s_date<check_time+1);

    if(isempty(pos))
        continue;
    end

    s_chla(s_chla<0) = NaN;

    figure(1);
    m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
    m_scatter(s_lon(pos),s_lat(pos),50,s_chla(pos),'filled');
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