clear all;close all;
addpath(path,'C:\Users\cheny\Desktop\matlab_tools')
load('OLD_ECOHAB_dat.mat');

filename = "C:\Users\cheny\Desktop\matlab_tools\m_map\private\gshhs_h.b";
S = gshhs(filename{1},[26.2 27.3],[-82.9 -81.9]);
clat = [S.Lat];
clon = [S.Lon];

[dat_year,~,~,~,~,~] = datevec(dat_date);
pos11 = find(dat_year==2008&~isnan(dat_syn));
pos21= find(dat_year==2009&~isnan(dat_syn));
pos31 = find(dat_year==2010&~isnan(dat_syn));

pos12 = find(dat_year==2008&~isnan(dat_kb));
pos22= find(dat_year==2009&~isnan(dat_kb));
pos32 = find(dat_year==2010&~isnan(dat_kb));

%%
figure('Units','pixels','Position',[100 100 1100 700]);
tiledlayout(2,3);
nexttile;
colormap("jet");
scatter(dat_lon(pos11),dat_lat(pos11),60,log10(dat_syn(pos11)*1000),'filled','MarkerEdgeColor',[0.7 0.7 0.7],'LineWidth',1.5);
t = colorbar;
ylabel(t,'log_1_0 Cells L^-^1');
hold on;
fill_coastline(clon,clat);
axis image
axis([-83.5 -81.2 25.4 28.2]);
set(gcf,'color','w')
caxis([6.5 8]);

nexttile;
colormap("jet");
scatter(dat_lon(pos21),dat_lat(pos21),60,log10(dat_syn(pos21)*1000),'filled','MarkerEdgeColor',[0.9 0.9 0.9],'LineWidth',1.5);
t = colorbar;
ylabel(t,'log_1_0 Cells L^-^1');
hold on;
fill_coastline(clon,clat);
axis image
axis([-83.5 -81.2 25.4 28.2]);
set(gcf,'color','w')
caxis([6.5 8]);

nexttile;
colormap("jet");
scatter(dat_lon(pos31),dat_lat(pos31),60,log10(dat_syn(pos31)*1000),'filled','MarkerEdgeColor','w','LineWidth',1.5);
t = colorbar;
ylabel(t,'log_1_0 Cells L^-^1');
hold on;
fill_coastline(clon,clat);
axis image
axis([-83.5 -81.2 25.4 28.2]);
set(gcf,'color','w')
caxis([6.5 8]);

nexttile;
colormap("jet");
scatter(dat_lon(pos12),dat_lat(pos12),60,log10(dat_kb(pos12)),'filled','MarkerEdgeColor',[0.7 0.7 0.7],'LineWidth',1.5);
t = colorbar;
ylabel(t,'log_1_0 Cells L^-^1');
hold on;
fill_coastline(clon,clat);
axis image
axis([-83.5 -81.2 25.4 28.2]);
set(gcf,'color','w')
caxis([3 6]);

nexttile;
colormap("jet");
scatter(dat_lon(pos22),dat_lat(pos22),60,log10(dat_kb(pos22)),'filled','MarkerEdgeColor',[0.9 0.9 0.9],'LineWidth',1.5);
t = colorbar;
ylabel(t,'log_1_0 Cells L^-^1');
hold on;
fill_coastline(clon,clat);
axis image
axis([-83.5 -81.2 25.4 28.2]);
set(gcf,'color','w')
caxis([3 6]);

nexttile;
colormap("jet");
scatter(dat_lon(pos32),dat_lat(pos32),60,log10(dat_kb(pos32)),'filled','MarkerEdgeColor','w','LineWidth',1.5);
t = colorbar;
ylabel(t,'log_1_0 Cells L^-^1');
hold on;
fill_coastline(clon,clat);
axis image
axis([-83.5 -81.2 25.4 28.2]);
set(gcf,'color','w')
caxis([3 6]);