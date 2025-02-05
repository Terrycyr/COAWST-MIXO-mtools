clear all;close all;
addpath(path,'../SUNRISE_SUNSET/');
addpath C:\Users\cheny\Desktop\EcoHAB\self_functions

grd = '../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd,'Cs_r'));

year = 2005;
out_date = datenum(year,1,1,0,0,0):3/24:datenum(year,12,31,24,0,0);
method = 2; %1: linear interpolate par; 2: use climatology

ndays= datenum(year+1,1,1,0,0,0)-datenum(year,1,1);

[r,c] = size(lon);

for i = 1:ndays
    dnum = i+datenum(year,1,1)-1;
    [rise(i),set(i),noon(i)] = sunrise(27.05,-82.75,0,0,datestr(dnum,'yyyy-mm-dd'));
    noon(i) = noon(i) - datenum(year,1,1);
    F(i) = set(i)-rise(i);
    rise(i) = rise(i) - datenum(year,1,1);
    set(i) = set(i) - datenum(year,1,1);
end

if(method==1)

    load(['./PAR_raw/PAR_raw_',num2str(year),'.mat']);
    tvec = datevec(time);
    par0 = par;
    clear par;

    for i = 1:ndays
        z = squeeze(par0(i,:,:));
        tmp = griddata(plon(~isnan(z)),plat(~isnan(z)),z(~isnan(z)),lon,lat,'linear');
        tmp(isnan(tmp)) = griddata(plon(~isnan(z)),plat(~isnan(z)),z(~isnan(z)),lon(isnan(tmp)),lat(isnan(tmp)),'nearest');
        par_m(:,:,i) = tmp;
        clear tmp
    end

    par = zeros(r,c,length(out_date));
    for i=1:length(out_date)
        tday(i) = out_date(i) - datenum(year,1,1);
        day = find((tday(i)-rise)>0&(tday(i)-set)<0);
        if(~isempty(day))
            tmp = par_m(:,:,day)/F(day)*pi/2* sin(3.14159*(tday(i)-rise(day))/F(day));
            tmp(mask==0) = 0.0;
            par(:,:,i) = tmp;
            par_avg(i) = mean(mean(tmp(mask==1)));
            par_var(i) = sqrt(var(tmp(mask==1),0,'all'));
            clear tmp
        else
            par(:,:,i) = 0.0;
            par_avg(i) = 0.0;
            par_var(i) = 0.0;
        end
    end

    figure;
    plot(tday,par_avg);
    hold on;
    for i=1:ndays
        tmp(i) = mean([mean(par_avg([1:8]+8*(i-1))),mean(par_avg([2:9]+8*(i-1)))]);
        tmp2(i) = mean([mean(par_var([1:8]+8*(i-1))),mean(par_var([2:9]+8*(i-1)))]);
    end
    hold on;
    errorbar(0.5:(ndays-0.5),tmp,tmp2,'LineWidth',2);
    hold on;
    alpha(0.7);
    xlim([0 366]);


elseif(method==2)

    load(['./PAR_raw_2003_2022.mat']);
    tvec = datevec(time);
    par0 = par;
    clear par;

    for i = 1:ndays
        tmp = datevec(i-1+datenum(year,1,1));
        mon_clm = tmp(2);
        day_clm = tmp(3);
        pos_clm = find(tvec(:,2)==mon_clm&tvec(:,3)==day_clm);
        par_clm(:,:,i) = squeeze(mean(par0(pos_clm,:,:),1,'omitnan'));

        z = squeeze(par_clm(:,:,i));
        tmp = griddata(plon(~isnan(z)),plat(~isnan(z)),z(~isnan(z)),lon,lat,'linear');
        tmp(isnan(tmp)) = griddata(plon(~isnan(z)),plat(~isnan(z)),z(~isnan(z)),lon(isnan(tmp)),lat(isnan(tmp)),'nearest');
        par_m(:,:,i) = tmp;
        clear tmp
    end

    par = zeros(r,c,length(out_date));
    for i=1:length(out_date)
        tday(i) = out_date(i) - datenum(year,1,1);
        day = find((tday(i)-rise)>0&(tday(i)-set)<0);
        if(~isempty(day))
            tmp = par_m(:,:,day)/F(day)*pi/2* sin(3.14159*(tday(i)-rise(day))/F(day));
            tmp(mask==0) = 0.0;
            par(:,:,i) = tmp;
            par_avg(i) = mean(mean(tmp(mask==1)));
            par_var(i) = sqrt(var(tmp(mask==1),0,'all'));
            clear tmp
        else
            par(:,:,i) = 0.0;
            par_avg(i) = 0.0;
            par_var(i) = 0.0;
        end
    end

    figure;
    plot(tday,par_avg);
    hold on;
    for i=1:ndays
        tmp(i) = mean([mean(par_avg([1:8]+8*(i-1))),mean(par_avg([2:9]+8*(i-1)))]);
        tmp2(i) = mean([mean(par_var([1:8]+8*(i-1))),mean(par_var([2:9]+8*(i-1)))]);
    end
    hold on;
    errorbar(0.5:(ndays-0.5),tmp,tmp2,'LineWidth',2);
    hold on;
    alpha(0.7);
    xlim([0 366]);

%     for i=1:2921
%         figure(2);
%         contourf(lon,lat,par(:,:,i));
%     end
end


PAR = par;
save(['PAR_',num2str(year),'.mat'],'tday','PAR');

