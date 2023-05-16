function [] = gen_bay_nutrients_ini(cruise_dir,init_fn,grd,bay_mask,year,include_cruise,dep_range,time_range,time_range2)

addpath(path,'../../matlab_tools');

% model grid
fn = grd;

% params.
N= length(ncread(fn,'Cs_r'));
Vtransform = ncread(fn,'Vtransform');
Vstretching = ncread(fn,'Vstretching');
THETA_S = ncread(fn,'theta_s');
THETA_B = ncread(fn,'theta_b');
TCLINE = ncread(fn,'Tcline');
hc = ncread(fn,'hc');

%DISSOLVED & PARTICULATE
OCDP_r = 0.84;
ONDP_r = 0.60;
OPDP_r = 0.22;

%LABILE & REFRACTORY
OCLR_r = 0.35;
ONLR_r = 0.35;
OPLR_r = 0.35;


lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
h = ncread(fn,'h');

%run([cruise_dir,'/read_cruise_data.m']);
[cruise_ini] = process_cruise_data(cruise_dir,grd,time_range,time_range2);


%no23
[~,~,raw] = xlsread(init_fn,1);

k = 0;
for i=2:size(raw,1)

    switch raw{i,9}
        case 'ug/l'
            k=k+1;
            tmp(k) = raw{i,8}/1000;
        case 'mg/l'
            k=k+1;
            tmp(k) = raw{i,8};
        otherwise
            k=k+1;
            tmp(k) = raw{i,8}/1000;
    end

    no23(k) = tmp(k);
    lon_no23(k) = raw{i,3};
    lat_no23(k) = raw{i,2};
end

%no23_bay = barnes(lon_no23,lat_no23,no23_m,lon,lat,0.1,0.1);
no23_bay0 = griddata(lon_no23(:),lat_no23(:),no23(:),lon,lat);
no23_bay0(bay_mask==0) = NaN;

for i=1:N
    if(include_cruise(1))
        tmp = cruise_ini.no23_int(:,:,i);
        tmp(~isnan(no23_bay0)) = no23_bay0(~isnan(no23_bay0));
        no23_bay(:,:,i) = tmp;
    else
        no23_bay(:,:,i) = no23_bay0;
    end
end

figure;
%contourf(lon,lat,no23_bay(:,:,1),'linestyle','none');
hp = pcolor(lon,lat,no23_bay(:,:,1));
set(hp,'linestyle','none');
hold on;
scatter(lon_no23(:),lat_no23(:),30,no23(:),'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
caxis([0 0.05]);
title('NO3');

%nh34
[~,~,raw] = xlsread(init_fn,2);
k = 0;
for i=2:size(raw,1)

    switch raw{i,9}
        case 'ug/l'
            k=k+1;
            tmp(k) = raw{i,8}/1000;
        case 'mg/l'
            k=k+1;
            tmp(k) = raw{i,8};
        otherwise
            k=k+1;
            tmp(k) = raw{i,8}/1000;
    end

    nh4(k) = tmp(k);
    lon_nh4(k) = raw{i,3};
    lat_nh4(k) = raw{i,2};
end

nh4_bay0 = griddata(lon_nh4,lat_nh4,nh4,lon,lat);
nh4_bay02 = nh4_bay0;
nh4_bay02(isnan(nh4_bay0)) = griddata(lon_nh4,lat_nh4,nh4,lon(isnan(nh4_bay0)),lat(isnan(nh4_bay0)),'nearest');
nh4_bay0(bay_mask==0) = NaN;
nh4_bay02(bay_mask==0) = NaN;

for i=1:N
    nh4_bay(:,:,i) = nh4_bay0;  
    nh4_bay2(:,:,i) = nh4_bay02; 
end

figure;
%contourf(lon,lat,nh4_bay(:,:,1),'linestyle','none');
%hold on;
hp = pcolor(lon,lat,nh4_bay(:,:,1));
set(hp,'linestyle','none');
hold on;
scatter(lon_nh4,lat_nh4,30,nh4,'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
colorbar;
title('NH4');

%po4
[dat,txt,raw] = xlsread(init_fn,3);
k = 0;
for i=2:size(raw,1)

    switch raw{i,9}
        case 'ug/l'
            k=k+1;
            tmp(k) = raw{i,8}/1000;
        case 'mg/l'
            k=k+1;
            tmp(k) = raw{i,8};
        otherwise
            k=k+1;
            tmp(k) = raw{i,8}/1000;
    end

    po4(k) = tmp(k);
    lon_po4(k) = raw{i,3};
    lat_po4(k) = raw{i,2};
end

po4_bay0 = griddata(lon_po4(:),lat_po4(:),po4(:),lon,lat);
po4_bay0(bay_mask==0) = NaN;

for i=1:N
    if(include_cruise(2))
        tmp = cruise_ini.po4_int(:,:,i);
        tmp(~isnan(po4_bay0)) = po4_bay0(~isnan(po4_bay0));
        po4_bay(:,:,i) = tmp;
    else
        po4_bay(:,:,i) = po4_bay0;  
    end
end

figure;
hp = pcolor(lon,lat,po4_bay(:,:,1));
set(hp,'linestyle','none');
%contourf(lon,lat,po4_bay(:,:,1),'linestyle','none');
hold on;
scatter(lon_po4(:),lat_po4(:),30,po4(:),'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
title('PO4');

%sit
[dat,txt,raw] = xlsread(init_fn,4);
k = 0;
for i=2:size(raw,1)

    switch raw{i,9}
        case 'ug/l'
            k=k+1;
            tmp(k) = raw{i,8}/1000;
        case 'mg/l'
            k=k+1;
            tmp(k) = raw{i,8};
        otherwise
            k=k+1;
            tmp(k) = raw{i,8}/1000;
    end

    sit(k) = tmp(k);
    lon_sit(k) = raw{i,3};
    lat_sit(k) = raw{i,2};
end

sit_bay0 = griddata(lon_sit(:),lat_sit(:),sit(:),lon,lat);
sit_bay0(bay_mask==0) = NaN;

for i=1:N
    if(include_cruise(3))
        tmp = cruise_ini.sit_int(:,:,i);
        tmp(~isnan(sit_bay0)) = sit_bay0(~isnan(sit_bay0));
        sit_bay(:,:,i) = tmp;
    else
        sit_bay(:,:,i) = sit_bay0;
    end
end

figure;
hp = pcolor(lon,lat,sit_bay(:,:,1));
set(hp,'linestyle','none');
%contourf(lon,lat,sit_bay(:,:,1),'linestyle','none');
hold on;
scatter(lon_sit(:),lat_sit(:),30,sit(:),'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
title('SIT');

%oc
[dat,txt,raw] = xlsread(init_fn,5);
k = 0;
for i=2:size(raw,1)

    switch raw{i,9}
        case 'ug/l'
            k=k+1;
            tmp(k) = raw{i,8}/1000;
        case 'mg/l'
            k=k+1;
            tmp(k) = raw{i,8};
        otherwise
            k=k+1;
            tmp(k) = raw{i,8}/1000;
    end

    oc(k) = tmp(k);
    lon_oc(k) = raw{i,3};
    lat_oc(k) = raw{i,2};
end

oc_bay0 = griddata(lon_oc,lat_oc,oc,lon,lat);
%!!!!!!!!!!!!!!!!!!!!!!!!!
oc_bay0(:,:) = 3;
%!!!!!!!!!!!!!!!!!!!!!!!!!
oc_bay0(bay_mask==0) = NaN;
%oc_dep = griddata(lon,lat,h,lon_oc,lat_oc);

for i=1:N
    oc_bay(:,:,i) = oc_bay0;  
end

figure;
hp = pcolor(lon,lat,oc_bay(:,:,1));
set(hp,'linestyle','none');
%contourf(lon,lat,oc_bay(:,:,1),'linestyle','none');
hold on;
scatter(lon_oc,lat_oc,30,oc,'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
title('OC');

%don
% [dat,txt,raw] = xlsread(init_fn,6);
% 
% k = 0;
% for i=2:size(raw,1)
% 
%     switch raw{i,9}
%         case 'ug/l'
%             k=k+1;
%             tmp(k) = raw{i,8}/1000;
%         case 'mg/l'
%             k=k+1;
%             tmp(k) = raw{i,8};
%         otherwise
%             k=k+1;
%             tmp(k) = raw{i,8}/1000;
%     end
% 
%     on(k) = tmp(k);
%     lon_on(k) = raw{i,3};
%     lat_on(k) = raw{i,2};
% end
% 
% don = on*ONDP_r;
% 
% don_bay0 = griddata(lon_on(:),lat_on(:),don(:),lon,lat);
% don_bay0(mask==0) = 0;
% 
% for i=1:N
%     if(include_cruise(4))
%         tmp = cruise_ini.don_int(:,:,i);
%         tmp(~isnan(don_bay0)) = don_bay0(~isnan(don_bay0));
%         don_bay(:,:,i) = tmp;
%     else
%         don_bay(:,:,i) = don_bay0;
%     end
% end

%tn
[~,~,raw] = xlsread(init_fn,7);
k = 0;
for i=2:size(raw,1)

    switch raw{i,9}
        case 'ug/l'
            k=k+1;
            tmp(k) = raw{i,8}/1000;
        case 'mg/l'
            k=k+1;
            tmp(k) = raw{i,8};
        otherwise
            k=k+1;
            tmp(k) = raw{i,8}/1000;
    end

    tn(k) = tmp(k);
    lon_tn(k) = raw{i,3};
    lat_tn(k) = raw{i,2};
end

tn_bay0 = griddata(lon_tn,lat_tn,tn,lon,lat);
tn_bay0(bay_mask==0) = NaN;

for i=1:N
    tn_bay(:,:,i) = tn_bay0; 
    on_bay(:,:,i) = tn_bay(:,:,i)-no23_bay(:,:,i)-nh4_bay2(:,:,i);
    don_bay(:,:,i) = on_bay(:,:,i)*ONDP_r;
end

figure;
contourf(lon,lat,don_bay(:,:,1),'linestyle','none');
colorbar;
%hp = pcolor(lon,lat,don_bay(:,:,1));
%set(hp,'linestyle','none');
colorbar;
title('DON');

figure;
contourf(lon,lat,tn_bay(:,:,1),'linestyle','none');
hold on;
scatter(lon_tn,lat_tn,30,tn,'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
%hp = pcolor(lon,lat,tn_bay);
%set(hp,'linestyle','none');
colorbar;
title('TN');

%sal
[dat,txt,raw] = xlsread(init_fn,8);
k = 0;
for i=2:size(raw,1)
    k=k+1;
    sal(k) = raw{i,8};
    lon_sal(k) = raw{i,3};
    lat_sal(k) = raw{i,2};
end
sal_bay0 = griddata(lon_sal,lat_sal,sal,lon,lat);
sal_bay02 = sal_bay0;
sal_bay02(isnan(sal_bay0)) = griddata(lon_sal,lat_sal,sal,lon(isnan(sal_bay0)),lat(isnan(sal_bay0)),'nearest');
sal_bay0(bay_mask==0) = NaN;
sal_bay02(bay_mask==0) = NaN;

for i=1:N
    sal_bay(:,:,i) = sal_bay02;  
end

figure;
hp = pcolor(lon,lat,sal_bay(:,:,1));
set(hp,'linestyle','none');
%contourf(lon,lat,sal_bay(:,:,1),'linestyle','none');
hold on;
scatter(lon_sal(:),lat_sal(:),30,sal(:),'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
caxis([20 35]);
title('SAL');

%CHLA
[dat,txt,raw] = xlsread(init_fn,9);
k = 0;
for i=2:size(raw,1)
    switch raw{i,9}
        case 'ug/l'
            k=k+1;
            tmp(k) = raw{i,8}/1000;
        case 'mg/l'
            k=k+1;
            tmp(k) = raw{i,8};
        otherwise
            k=k+1;
            tmp(k) = raw{i,8}/1000;
    end

    chla(k) = tmp(k);
    lon_chla(k) = raw{i,3};
    lat_chla(k) = raw{i,2};
end
chla_bay0 = griddata(lon_chla,lat_chla,chla,lon,lat);
chla_bay0(bay_mask==0) = NaN;

for i=1:N
    chla_bay(:,:,i) = chla_bay0;  
end

figure;
hp = pcolor(lon,lat,chla_bay(:,:,1));
set(hp,'linestyle','none');
%contourf(lon,lat,chla_bay(:,:,1),'linestyle','none');
hold on;
scatter(lon_chla(:),lat_chla(:),30,chla(:),'filled','MarkerEdgeColor','r','LineWidth',1);
colorbar;
title('CHLA');

save(strcat('bay_ini_',num2str(year),'.mat'),'no23_bay','nh4_bay','po4_bay'...
    ,'sit_bay','oc_bay','don_bay','tn_bay','sal_bay','chla_bay');
end









