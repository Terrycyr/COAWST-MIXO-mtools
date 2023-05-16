function [cruise_ini] = process_cruise_data(cruise_dir,grd,time_range,time_range2)
%process_cruise_data Summary of this function goes here
%   Detailed explanation goes here


load([cruise_dir, '/Inorganicnutrients.mat']);
load([cruise_dir, '/Organicnutrients.mat']');

lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
Cs_r = ncread(grd,'Cs_r');
dep = ncread(grd,'h');
[r,c] = size(mask);
NZ = length(Cs_r);

pos = IN_num(:,1)>=time_range(1)&IN_num(:,1)<=time_range(2);
pos2 = DON_num(:,1)>=time_range2(1)&DON_num(:,1)<=time_range2(2);

lat_cin = IN_num(pos,2);
lon_cin = IN_num(pos,3);
station_cin = IN_out(pos,2);
dep_cin = IN_num(pos,4);
sit_cin = IN_num(pos,5)*28/1000;
no23_cin = IN_num(pos,6)*14/1000;
po4_cin = IN_num(pos,7)*31/1000;

lat_con = DON_num(pos2,2);
lon_con = DON_num(pos2,3);
station_con = DON_out(pos2,2);
dep_con = DON_num(pos2,4);
don_con = DON_num(pos2,5)*14/1000;

k = 0;
for i= 1:length(dep_cin)
    if(isnan(dep_cin(i)))
        k=k+1;
        del_pos(k) = i;
    end
end

k = 0;
for i= 1:length(dep_con)
    if(isnan(dep_con(i)))
        k=k+1;
        del_pos2(k) = i;
    end
end

if(exist('del_pos','var'))
lat_cin(del_pos) = [];
lon_cin(del_pos) = [];
station_cin(del_pos) = [];
sit_cin(del_pos) = [];
no23_cin(del_pos) = [];
po4_cin(del_pos) = [];
dep_cin(del_pos) = [];
end

if(exist('del_pos2','var'))
lat_con(del_pos2) = [];
lon_con(del_pos2) = [];
station_con(del_pos2) = [];
sit_con(del_pos2) = [];
dep_con(del_pos2) = [];
end

cate = double(categorical(station_cin));
cate2 = double(categorical(station_con));
k1=0;
k2=0;
k3=0;
k4=0;
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
        x2 = dep(s_i,s_j)*(linspace(1/42,41/42,21));
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
        x2 = dep(s_i,s_j)*(linspace(1/42,41/42,21));
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
        x2 = dep(s_i,s_j)*(linspace(1/42,41/42,21));
        if(length(y(~isnan(y)))>1)
            tmp = interp1(x(~isnan(y)),y(~isnan(y)),x2);
            tmp(isnan(tmp)) = interp1(x(~isnan(y)),y(~isnan(y)),x2(isnan(tmp)),'nearest','extrap');
            sit_out{4}(k3,1:NZ) = tmp;
        else
            sit_out{4}(k3,1:NZ) = y(~isnan(y));
        end
    end
end

clear c_lon c_lat c_dep c_station c_pos
for i=1:max(cate2)
    c_pos = find(cate2==i);
    c_lon(i) = lon_con(c_pos(1));
    c_lat(i) = lat_con(c_pos(1));
    c_dep{i} = dep_con(c_pos);
    c_don{i} = don_con(c_pos);
    c_station{i} = station_cin(c_pos);

    dep_cate = double(categorical(c_dep{i}));
    for j=1:max(dep_cate)
        dep_pos = find(dep_cate==j);
        cdep{i}(j) = median(c_dep{i}(dep_pos),'omitnan');
        cdon{i}(j) = median(c_don{i}(dep_pos),'omitnan');
    end

    if(sum(~isnan(cdon{i}))&&sum(cdon{i}>0))
        k4=k4+1;
        don_out{1}(k4) = c_lon(i);
        don_out{2}(k4) = c_lat(i);
        don_out{3}{k4} = c_station{i}{1};
        [s_i,s_j] = find_ij(c_lon(i),c_lat(i),lon,lat,mask);
        x = cdep{i};
        y = cdon{i};
        x2 = dep(s_i,s_j)*(linspace(1/42,41/42,21));
        if(length(y(~isnan(y)))>1)
            tmp = interp1(x(~isnan(y)),y(~isnan(y)),x2);
            tmp(isnan(tmp)) = interp1(x(~isnan(y)),y(~isnan(y)),x2(isnan(tmp)),'nearest','extrap');
            don_out{4}(k4,1:NZ) = tmp;
        else
            don_out{4}(k4,1:NZ) = y(~isnan(y));
        end
    end
end

for i=1:NZ
    try
        no23_int(:,:,i) = griddata(no23_out{1},no23_out{2},no23_out{4}(:,i),lon,lat,'natural');
    catch
        no23_int(:,:,i) = zeros(size(lon));
    end
    try
        po4_int(:,:,i) = griddata(po4_out{1},po4_out{2},po4_out{4}(:,i),lon,lat,'natural');
    catch
        po4_int(:,:,i) = zeros(size(lon));
    end
    try
        sit_int(:,:,i) = griddata(sit_out{1},sit_out{2},sit_out{4}(:,i),lon,lat,'natural');
    catch
        sit_int(:,:,i) = zeros(size(lon));
    end
    try
        don_int(:,:,i) = griddata(don_out{1},don_out{2},don_out{4}(:,i),lon,lat,'natural');
    catch
        don_int(:,:,i) = zeros(size(lon));
    end
end

cruise_ini.no23_int = no23_int;
cruise_ini.po4_int = po4_int;
cruise_ini.sit_int = sit_int;
cruise_ini.don_int = don_int;

figure('Units','pixels','Position',[100 100 1000 800]);
colormap(jet);

tiledlayout(4,3,'TileSpacing', 'compact');

%no3
if(exist('no23_out','var'))
nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,no23_int(:,:,21),'linestyle','none');
hold on;
m_scatter(no23_out{1},no23_out{2},40,no23_out{4}(:,21),'filled','MarkerEdgeColor','w');
hold on;
m_text(no23_out{1}+0.08,no23_out{2},no23_out{3},'color','r','fontsize',7);
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('NO23');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,no23_int(:,:,10),'linestyle','none');
hold on;
m_scatter(no23_out{1},no23_out{2},40,no23_out{4}(:,10),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('NO23');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,no23_int(:,:,1),'linestyle','none');
hold on;
m_scatter(no23_out{1},no23_out{2},40,no23_out{4}(:,1),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('NO23');
end

%po4
if(exist('po4_out','var'))
nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,po4_int(:,:,21),'linestyle','none');
hold on;
m_scatter(po4_out{1},po4_out{2},40,po4_out{4}(:,21),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('PO4');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,po4_int(:,:,10),'linestyle','none');
hold on;
m_scatter(po4_out{1},po4_out{2},40,po4_out{4}(:,10),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('PO4');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,po4_int(:,:,1),'linestyle','none');
hold on;
m_scatter(po4_out{1},po4_out{2},40,po4_out{4}(:,1),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('PO4');
end

%sit
if(exist('sit_out','var'))
nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,sit_int(:,:,21),'linestyle','none');
hold on;
m_scatter(sit_out{1},sit_out{2},40,sit_out{4}(:,21),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('SIT');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,sit_int(:,:,10),'linestyle','none');
hold on;
m_scatter(sit_out{1},sit_out{2},40,sit_out{4}(:,10),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('SIT');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,sit_int(:,:,1),'linestyle','none');
hold on;
m_scatter(sit_out{1},sit_out{2},40,sit_out{4}(:,1),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('SIT');
end

%don
if(exist("don_out",'var'))
nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,don_int(:,:,21),'linestyle','none');
hold on;
m_scatter(don_out{1},don_out{2},40,don_out{4}(:,21),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('DON');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,don_int(:,:,10),'linestyle','none');
hold on;
m_scatter(don_out{1},don_out{2},40,don_out{4}(:,10),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('DON');

nexttile;
m_proj('miller','lat',[25 28.5],'long',[-85 -81]);
hold on;
m_gshhs_f('patch',[.8 .8 .8]);
m_grid('linestyle','none','linewidth',1,'fontsize',11);
m_contourf(lon,lat,don_int(:,:,1),'linestyle','none');
hold on;
m_scatter(don_out{1},don_out{2},40,don_out{4}(:,1),'filled','MarkerEdgeColor','w');
%caxis([0 75]);
colorbar;
set(gcf,'color','w');
title('DON');
end

end