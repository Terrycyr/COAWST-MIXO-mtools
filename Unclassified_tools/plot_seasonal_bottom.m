clear all; close all;

for season = 1:4

cruise_dir = './';

load([cruise_dir, '/Inorganicnutrients.mat']);
load([cruise_dir, '/Organicnutrients.mat']');
grd = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';

lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
dep = ncread(grd,'h');
[r,c] = size(mask);
NZ = 21;


    season
    switch season
        case 1
            time_range = [datenum(2000,12,1),datenum(2001,3,1)];
            pos = IN_num(:,1)>time_range(1)&IN_num(:,1)<=time_range(2);
        case 2
            time_range = [datenum(2000,3,1),datenum(2001,6,1)];
            pos = IN_num(:,1)>time_range(1)&IN_num(:,1)<=time_range(2);
        case 3
            time_range = [datenum(2000,6,1),datenum(2001,8,1)];
            pos = IN_num(:,1)>time_range(1)&IN_num(:,1)<=time_range(2);
        case 4
            time_range = [datenum(2000,9,1),datenum(2001,12,1)];
            pos = IN_num(:,1)>time_range(1)&IN_num(:,1)<=time_range(2);
    end

%no23
lat_cin = IN_num(pos,2);
lon_cin = IN_num(pos,3);
station_cin = IN_out(pos,2);
dep_cin = IN_num(pos,4);
sit_cin = IN_num(pos,5)*28/1000;
no23_cin = IN_num(pos,6)*14/1000;
po4_cin = IN_num(pos,7)*31/1000;

k = 0;
for i= 1:length(dep_cin)
    if(isnan(dep_cin(i)))
        k=k+1;
        del_pos(k) = i;
    end
end


try
    lat_cin(del_pos) = [];
    lon_cin(del_pos) = [];
    station_cin(del_pos) = [];
    sit_cin(del_pos) = [];
    no23_cin(del_pos) = [];
    po4_cin(del_pos) = [];
    dep_cin(del_pos) = [];
catch
end


cate = double(categorical(station_cin));
k1=0;
k2=0;
k3=0;
for i=1:max(cate)
    c_pos = find(cate==i);
    c_lon(i) = lon_cin(c_pos(1));
    c_lat(i) = lat_cin(c_pos(1));
    c_dep{i} = dep_cin(c_pos);
    c_sit{i} = sit_cin(c_pos);
    c_no23{i} = no23_cin(c_pos);
    c_po4{i} = po4_cin(c_pos);
    c_station{i} = station_cin(c_pos);

    dep_cate = double(categorical(c_dep{i}));
    for j=1:max(dep_cate)
        dep_pos = find(dep_cate==j);
        cdep{i}(j) = median(c_dep{i}(dep_pos),'omitnan');
        cno23{i}(j) = median(c_no23{i}(dep_pos),'omitnan');
        cpo4{i}(j) = median(c_po4{i}(dep_pos),'omitnan');
        csit{i}(j) = median(c_sit{i}(dep_pos),'omitnan');
    end

    if(sum(~isnan(cno23{i}))&&sum(cno23{i}>0))
        k1=k1+1;
        no23_out{1}(k1) = c_lon(i);
        no23_out{2}(k1) = c_lat(i);
        no23_out{3}{k1} = c_station{i}{1};
        [s_i,s_j] = find_ij(c_lon(i),c_lat(i),lon,lat,mask);
        x = cdep{i};
        y = cno23{i};
        x2 = dep(s_i,s_j)*(linspace(41/42,1/42,21));
        if(length(y(~isnan(y)))>1)
            tmp = interp1(x(~isnan(y)),y(~isnan(y)),x2);
            tmp(isnan(tmp)) = interp1(x(~isnan(y)),y(~isnan(y)),x2(isnan(tmp)),'nearest','extrap');
            no23_out{4}(k1,1:NZ) = tmp;
        else
            no23_out{4}(k1,1:NZ) = y(~isnan(y));
        end
    end

    if(sum(~isnan(cpo4{i}))&&sum(cpo4{i}>0))
        k2=k2+1;
        po4_out{1}(k2) = c_lon(i);
        po4_out{2}(k2) = c_lat(i);
        po4_out{3}{k2} = c_station{i}{1};
        [s_i,s_j] = find_ij(c_lon(i),c_lat(i),lon,lat,mask);
        x = cdep{i};
        y = cpo4{i};
        x2 = dep(s_i,s_j)*(linspace(41/42,1/42,21));
        if(length(y(~isnan(y)))>1)
            tmp = interp1(x(~isnan(y)),y(~isnan(y)),x2);
            tmp(isnan(tmp)) = interp1(x(~isnan(y)),y(~isnan(y)),x2(isnan(tmp)),'nearest','extrap');
            po4_out{4}(k2,1:NZ) = tmp;
        else
            po4_out{4}(k2,1:NZ) = y(~isnan(y));
        end
    end

    if(sum(~isnan(csit{i}))&&sum(csit{i}>0))
        k3=k3+1;
        sit_out{1}(k3) = c_lon(i);
        sit_out{2}(k3) = c_lat(i);
        sit_out{3}{k3} = c_station{i}{1};
        [s_i,s_j] = find_ij(c_lon(i),c_lat(i),lon,lat,mask);
        x = cdep{i};
        y = csit{i};
        x2 = dep(s_i,s_j)*(linspace(41/42,1/42,21));
        if(length(y(~isnan(y)))>1)
            tmp = interp1(x(~isnan(y)),y(~isnan(y)),x2);
            tmp(isnan(tmp)) = interp1(x(~isnan(y)),y(~isnan(y)),x2(isnan(tmp)),'nearest','extrap');
            sit_out{4}(k3,1:NZ) = tmp;
        else
            sit_out{4}(k3,1:NZ) = y(~isnan(y));
        end
    end
end

for i=1:NZ
    no23_int(:,:,i) = griddata(no23_out{1},no23_out{2},no23_out{4}(:,i),lon,lat,'natural');
    po4_int(:,:,i) = griddata(po4_out{1},po4_out{2},po4_out{4}(:,i),lon,lat,'natural');
    sit_int(:,:,i) = griddata(sit_out{1},sit_out{2},sit_out{4}(:,i),lon,lat,'natural');
end

if(season==1)
figure('Units','pixels','Position',[100 100 1000 800]);
colormap(jet);
tiledlayout(4,3,'TileSpacing', 'compact');
end

%no3
nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
hp=m_pcolor(lon,lat,no23_int(:,:,1)/14*1000);
set(hp,'linestyle','none');
hold on;
[C,h] = m_contour(lon,lat,dep,[0 20 40 60 100],'color',[0.6 0.6 0.6],'linestyle','--');
clabel(C,h);
%m_scatter(no23_out{1},no23_out{2},40,no23_out{4}(:,1),'filled','MarkerEdgeColor','w');
caxis([0 10]);
colorbar;
set(gcf,'color','w');

%po4

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
hp = m_pcolor(lon,lat,po4_int(:,:,1)/31*1000);
set(hp,'linestyle','none')
hold on;
[C,h] = m_contour(lon,lat,dep,[0 20 40 60 100],'color',[0.6 0.6 0.6],'linestyle','--');
clabel(C,h);
%m_scatter(po4_out{1},po4_out{2},40,po4_out{4}(:,1),'filled','MarkerEdgeColor','w');
caxis([0 1]);
colorbar;
set(gcf,'color','w');

%sit
nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
hp = m_pcolor(lon,lat,sit_int(:,:,1)/28*1000);
set(hp,'linestyle','none');
hold on;
[C,h] = m_contour(lon,lat,dep,[0 20 40 60 100],'color',[0.6 0.6 0.6],'linestyle','--');
clabel(C,h);
%m_scatter(sit_out{1},sit_out{2},40,sit_out{4}(:,1),'filled','MarkerEdgeColor','w');
caxis([0 8]);
colorbar;
set(gcf,'color','w');

clear all
end