clear all;close all;
addpath(path,'../SUNRISE_SUNSET/');
addpath C:\Users\cheny\Desktop\EcoHAB\self_functions

grd = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd,'Cs_r'));

year = 2021;
out_date = datenum(year,1,1,0,0,0):3/24:datenum(year,12,31,24,0,0);
method = 3; %1: linear interpolate par; 2: swrad linear fitting to climatology 3: use climatology

ndays= datenum(year+1,1,1,0,0,0)-datenum(year,1,1);

fn1 = ['../NC_file_generation/WFS_',num2str(year),'_swrad1.nc'];

load('C:\Users\cheny\Desktop\EcoHAB\PAR\erdMGpar01day_c46c_e27d_dbdf.mat')

time = erdMGpar01day.time;
tvec = datevec(time/3600/24+datenum(1970,1,1));

par_lon = mean(erdMGpar01day.longitude);
par_lat = mean(erdMGpar01day.latitude);
par0 = mean(mean(double(erdMGpar01day.par),4),3);

for i = 1:ndays
    dnum = i+datenum(year,1,1)-1;
    [rise(i),set(i),noon(i)] = sunrise(27.05,-82.75,0,0,datestr(dnum,'yyyy-mm-dd'));
    noon(i) = noon(i) - datenum(year,1,1);
    F(i) = set(i)-rise(i);
    rise(i) = rise(i) - datenum(year,1,1);
    set(i) = set(i) - datenum(year,1,1);
end

if(method==1)

    x0 = time/3600/24+datenum(1970,1,1);
    y0 = par0*1e6/24/3600 ;%Einsteins m-2 d-1 -> µmoles m-2 s-1

    par_m = interp1(x0(~isnan(y0)),y0(~isnan(y0)),[1:ndays]+datenum(year,1,0,12,0,0));

    for i=1:length(out_date)
        tday(i) = out_date(i) - datenum(year,1,1);
        day = find((tday(i)-rise)>0&(tday(i)-set)<0);
        if(~isempty(day))
            par(i) = par_m(day)/F(day)*pi/2* sin(3.14159*(tday(i)-rise(day))/F(day));
        else
            par(i) = 0.0;
        end
    end

    figure;
    plot(tday,par);
    hold on;
    for i=1:ndays
        tmp(i) = mean([mean(par([1:8]+8*(i-1))),mean(par([2:9]+8*(i-1)))]);
    end
    hold on;
    scatter(x0-datenum(year,1,1),y0,'filled')
    hold on;
    plot(0.5:(ndays-0.5),tmp,'LineWidth',2)
    alpha(0.7);
    xlim([0 366]);

elseif(method==2)

    for i = 1:ndays
        tmp = datevec(i-1+datenum(year,1,1));
        mon_clm = tmp(2);
        day_clm = tmp(3);
        pos_clm = find(tvec(:,2)==mon_clm&tvec(:,3)==day_clm);
        date_clm(i) = mean(time(pos_clm)/3600/24+datenum(1970,1,1));
        par_clm(i) = mean(par0(pos_clm),'omitnan')*1e6/24/3600; %Einsteins m-2 d-1 -> µmoles m-2 s-1
    end

    [i,j] = find_ij(par_lon,par_lat,lon,lat,mask);

    swrad = squeeze(double(ncread(fn1,'swrad',[i,j,1],[1,1,inf])));

    for i=1:ndays
        swrad_m(i) = mean([mean(swrad([1:4]+4*(i-1))),mean(swrad([2:5]+4*(i-1)))]);
    end

    %Linear Fit swrad_m v.s.par_m
    [xData, yData] = prepareCurveData( swrad_m, par_clm );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'LAR';

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    %par_m = fitresult(swrad_m);
    par_m = swrad_m*1.903+89.12;

    for i=1:length(out_date)
        tday(i) = out_date(i) - datenum(year,1,1);
        day = find((tday(i)-rise)>0&(tday(i)-set)<0);
        if(~isempty(day))
            par(i) = par_m(day)/F(day)*pi/2* sin(3.14159*(tday(i)-rise(day))/F(day));
        else
            par(i) = 0.0;
        end
    end


    figure;
    plot(tday,par);
    hold on;
    for i=1:ndays
        tmp(i) = mean([mean(par([1:8]+8*(i-1))),mean(par([2:9]+8*(i-1)))]);
    end
    hold on;
    scatter(0.5:(ndays-0.5),par_m,'g','filled');
    hold on;
    plot(0.5:(ndays-0.5),tmp,'r','LineWidth',2)
    legend('MODEL','FITTED MEAN','MODEL CALCULATED MEAN ');
    alpha(0.7);
    xlim([100 120]);

    figure;
    plot(tday,par);
    hold on;
    scatter(0.5:(ndays-0.5),par_m,'g','filled');
    hold on;
    scatter(0.5:(ndays-0.5),par_clm(:),'k','filled');
    alpha(0.3);
    hold on;
    plot(0.5:(ndays-0.5),tmp,'r','LineWidth',2);
    legend('MODEL','FITTED MEAN','CLM','MODEL CALCULATED MEAN ');

elseif(method==3)

    for i = 1:ndays
        tmp = datevec(i-1+datenum(year,1,1));
        mon_clm = tmp(2);
        day_clm = tmp(3);
        pos_clm = find(tvec(:,2)==mon_clm&tvec(:,3)==day_clm);
        date_clm(i) = mean(time(pos_clm)/3600/24+datenum(1970,1,1));
        par_clm(i) = mean(par0(pos_clm),'omitnan')*1e6/24/3600; %Einsteins m-2 d-1 -> µmoles m-2 s-1
    end

    for i=1:length(out_date)
        tday(i) = out_date(i) - datenum(year,1,1);
        day = find((tday(i)-rise)>0&(tday(i)-set)<0);
        if(~isempty(day))
            par(i) = par_clm(day)/F(day)*pi/2* sin(3.14159*(tday(i)-rise(day))/F(day));
        else
            par(i) = 0.0;
        end
    end

    figure;
    plot(tday,par);
    hold on;
    scatter(0.5:(ndays-0.5),par_clm(:),'k','filled');
    alpha(0.3);
    hold on;
    for i=1:ndays
        tmp(i) = mean([mean(par([1:8]+8*(i-1))),mean(par([2:9]+8*(i-1)))]);
    end
    plot(0.5:(ndays-0.5),tmp,'r','LineWidth',2);
    legend('MODEL','CLM','MODEL CALCULATED MEAN ');
end


ITOTSF = par;
save(['PAR_',num2str(year),'.mat'],'tday','ITOTSF');

