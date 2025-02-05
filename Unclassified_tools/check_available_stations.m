clear all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

[dat,txt,raw] = xlsread('Check.xlsx');

for i=2:123
    loc(i-1,1) = raw{i,7};
    loc(i-1,2) = raw{i,8};
    loc_name{i-1} = raw{i,6};
    loc_time(i-1) = xlstime2date(raw{i,10});
end

loc_tvec = datevec(loc_time);
loc_mon = loc_tvec(:,2);
%%
figure(1);
m_proj('Mercator','lat',[25.8 27.3],'long',[-82.9 -81.6]);
%m_scatter(-82.059632,26.950358,'r','filled');
hold on;
for i=1:size(loc,1)
    i
    m_scatter(loc(i,2),loc(i,1),40,loc_mon(i),'filled');
    hold on;
    m_text(loc(i,2),loc(i,1),[datestr(loc_time(i),'mm-yy') loc_name{i}]);
end

m_gshhs_h('patch',[.6 .6 .6]);
hold on;
m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
colorbar;
set(gcf,'color', [1 1 1]);
xlabel('Longitude','fontsize',9);
ylabel('Latitude','fontsize',9);
hold off;

figure(2);
for i=1:size(loc,1)
scatter(loc(i,2),loc(i,1),40,loc_mon(i),'filled');
hold on;
text(loc(i,2),loc(i,1),[datestr(loc_time(i),'mm-yy') loc_name{i}]);
end
