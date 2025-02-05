clear all;close all;
addpath(path,'../../matlab_tools/m_map/');
addpath(path,'../self_functions/');
[~,~,raw] = xlsread('TB2021.xlsx',1);

st_date = xlstime2date(cell2mat(raw(2:end,13)));
st_tvec = datevec(st_date);
st_lat = cell2mat(raw(2:end,10));
st_lon = cell2mat(raw(2:end,11));
st_bdo = cell2mat(raw(2:end,41));

figure('Units','pixels','Position',[100 100 800 600]);
colormap(jet);
tiledlayout(3,3,"Padding","compact","TileSpacing","compact");
for i=4:12
    pos = find(st_tvec(:,2)==i);
    x = st_lon(pos);
    y = st_lat(pos);
    v = st_bdo(pos);
    nexttile;


    m_proj('Mercator','lat',[27.4 28.2],'long',[-82.9 -82.2]);
    m_scatter(x,y,40,v,'filled');
    m_gshhs_h('patch',[.8 .8 .8],'linestyle','none');
    if(i~=10)
        m_grid('linestyle','none','xtick',[-82.8:0.2:-82.2],'ytick',[-27.5:0.2:28.1],'fontsize',11,'tickstyle','dd','XTickLabel',[],'YTickLabel',[]);
    else
        m_grid('linestyle','none','xtick',[-82.8:0.2:-82.2],'ytick',[-27.5:0.2:28.1],'fontsize',11,'tickstyle','dd');
    end
    
    title(datestr(datenum(2021,i,1),'mmm'));
    set(gcf,'color','w');

    caxis([2 8]);
    if(i==12)
        colorbar;
    end
end