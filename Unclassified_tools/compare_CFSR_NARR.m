clear all; close all;

addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\Atmosphere_forcing_preprocessing');

%37.032 N 76.083 W CHBV2 - 8638901
%36.609 N 74.842 W 44014

station_lon = [-76.083 -74.842];
station_lat = [37.032 36.609];


dat1 = dlmread("chbv2h2022.txt",' ',2,0);
dat1(dat1==99) =  NaN;
year{1} = dat1(:,1);
mon{1} = dat1(:,2);
day{1} = dat1(:,3);
hour{1} = dat1(:,4);
min{1} = dat1(:,5);
wind_date{1} = datenum(year{1},mon{1},day{1},hour{1},min{1},0);
tmp = dat1(:,8);
tmp(tmp==0) = NaN;
wind_speed{1} = tmp;
tmp = dat1(:,6);
tmp(tmp==0) = NaN;
tmp(tmp>360) = NaN;
wind_dir{1} = tmp;


dat2 = dlmread("44014h2022.txt",' ',2,0);
dat2(dat2==99) =  NaN;
year{2} = dat2(:,1);
mon{2} = dat2(:,2);
day{2} = dat2(:,3);
hour{2} = dat2(:,4);
min{2} = dat2(:,5);
wind_date{2} = datenum(year{2},mon{2},day{2},hour{2},min{2},0);
tmp = dat2(:,8);
tmp(tmp==0) = NaN;
wind_speed{2} = tmp;
tmp = dat2(:,6);
tmp(tmp==0) = NaN;
tmp(tmp>360) = NaN;
wind_dir{2} = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify Variables
cfsr_dir = {'./','./'}; %CFSR directory for each variables
cfsr_var = {'U_GRD_L103','V_GRD_L103'}; %CFSR variable names
cfsr_deltat = 6; %time interval of the CFSR data, unit: hour
out_var = {'u_out','v_out'}; %name of the output variables
origin_date = datenum(2022,1,1,6,0,0):cfsr_deltat/24:datenum(2022,12,31,24,0,0); %origin_date should Cover the Initial time
out_date = origin_date;
read_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------NO NEED TO CHANGE BELOW--------------------------
year = datevec(out_date(1));
year = year(1);
outname = strcat('wind_',num2str(year),'.mat');
outname2 = strcat('wind_',num2str(year),'_NARR.mat');
lon = station_lon;
lat = station_lat;

if(read_flag==1)
    %Read data
    read_dat = cell(1,length(out_var));
    lon_cfsr = cell(1,length(out_var));
    lat_cfsr = cell(1,length(out_var));
    for time = origin_date(1):cfsr_deltat/24:origin_date(end)
        for var_n = 1:length(out_var)
            [cfsr_fname,time_id] = filename_tid(cfsr_dir{var_n},time,cfsr_var{var_n});
            tmp = ncread(cfsr_fname,cfsr_var{var_n},[1,1,time_id],[Inf Inf 1]);
            if(time == origin_date(1))
                x = double(ncread(cfsr_fname,'lon'));
                x(x>180) = x(x>180)-360;
                y = double(ncread(cfsr_fname,'lat'));
                [lat_cfsr{var_n},lon_cfsr{var_n}] = meshgrid(y,x);
                read_dat{var_n} = tmp;
            else
                tsize_tmp = size(tmp,3);
                tsize_dat = size(read_dat{var_n},3);
                read_dat{var_n}(:,:,[1:tsize_tmp]+tsize_dat) = tmp;
            end

        end
    end

    %Pre-processing for Saving
    tinterp_dat = cell(1,length(out_var));
    sinterp_dat = cell(1,length(out_var));
    for var_n = 1:length(out_var)
        %Temporal interpolation
        tinterp_dat{var_n} = tv3dinterpt(read_dat{var_n}, origin_date, out_date, 'linear');
        %Spatial interpolation
        input_x = lon_cfsr{var_n};
        input_y = lat_cfsr{var_n};
        out_x = lon;
        out_y = lat;
        sinterp_dat{var_n} = tv3dinterps(tinterp_dat{var_n}, input_x, input_y, out_x, out_y,'linear');
    end

    %Saving
    for var_n = 1:length(out_var)
        eval(strcat(out_var{var_n},'= sinterp_dat{var_n};'));
    end

    if(~exist(outname,'file'))
        save(outname,"out_date",'-v7.3')
    end
    for var_n = 1:length(out_var)
        eval(strcat('save(outname,','out_var{var_n},''-append'');'));
    end

%NARR DATA
    
    fn1 = './NARR/uwnd.10m.2022.nc';
    fn2 = './NARR/vwnd.10m.2022.nc';
    raw_time = ncread(fn1,'time')/24+datenum(1800,1,1);
    raw_lon = double(ncread(fn1,'lon'));
    raw_lat = double(ncread(fn1,'lat'));
    for st=1:length(station_lat)
        [st_i(st),st_j(st)] = find_ij(station_lon(st),station_lat(st),raw_lon,raw_lat,ones(size(raw_lon)));
        raw_u0(:,:,:,st) = double(squeeze(ncread(fn1,'uwnd',[st_i(st)-4 st_j(st)-4 1],[9 9 Inf])));
        raw_v0(:,:,:,st) = double(squeeze(ncread(fn2,'vwnd',[st_i(st)-4 st_j(st)-4 1],[9 9 Inf])));
        raw_lon2{st} = raw_lon(st_i(st)-4:st_i(st)+4,st_j(st)-4:st_j(st)+4);
        raw_lat2{st} = raw_lat(st_i(st)-4:st_i(st)+4,st_j(st)-4:st_j(st)+4);
        for i=1:length(raw_time)
            i
            raw_u(i,st) =  griddata(raw_lon2{st},raw_lat2{st},raw_u0(:,:,i,st),station_lon(st),station_lat(st));
            raw_v(i,st) =  griddata(raw_lon2{st},raw_lat2{st},raw_v0(:,:,i,st),station_lon(st),station_lat(st));
        end
    end
    
    narr_t = raw_time;
    narr_u = raw_u;
    narr_v = raw_v;

    save(outname2,'narr_v','narr_u','narr_t');
else
    load(outname);
    load(outname2);
end





cfsr_wind_speed = sqrt(squeeze(u_out.^2)+squeeze(v_out.^2));
cfsr_dir = 270 - (180/pi)*atan2(squeeze(v_out),squeeze(u_out));
for i = 1:size(cfsr_dir(1,:),2)
    tmp = cfsr_dir(1,:);
    tmp(tmp>360) = tmp(tmp>360)-360;
    cfsr_dir(1,:) = tmp;
end
for i = 1:size(cfsr_dir(2,:),2)
    tmp = cfsr_dir(2,:);
    tmp(tmp>360) = tmp(tmp>360)-360;
    cfsr_dir(2,:) = tmp;
end

narr_wind_speed0 = sqrt(squeeze(narr_u.^2)+squeeze(narr_v.^2));
narr_wind_speed(:,1) = interp1(narr_t,narr_wind_speed0(:,1),out_date);
narr_wind_speed(:,2) = interp1(narr_t,narr_wind_speed0(:,2),out_date);
narr_dir0 = 270 - (180/pi)*atan2(squeeze(narr_v),squeeze(narr_u));
for i = 1:size(narr_dir0(:,1))
    tmp = narr_dir0(:,1);
    tmp(tmp>360) = tmp(tmp>360)-360;
    narr_dir0(:,1) = tmp;
end
for i = 1:size(narr_dir0(:,2))
    tmp = narr_dir0(:,2);
    tmp(tmp>360) = tmp(tmp>360)-360;
    narr_dir0(:,2) = tmp;
end
narr_dir(:,1) = interp1(narr_t,narr_dir0(:,1),out_date);
narr_dir(:,2) = interp1(narr_t,narr_dir0(:,2),out_date);

buoy_wind(:,1) = interp1(wind_date{1},wind_speed{1},out_date);
buoy_wind(:,2) = interp1(wind_date{2},wind_speed{2},out_date);
buoy_dir(:,1) = interp1(wind_date{1},wind_dir{1},out_date);
buoy_dir(:,2) = interp1(wind_date{2},wind_dir{2},out_date);

figure('Units','pixels','Position',[100 100 800 500])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');


nexttile;
plot(out_date,buoy_wind(:,1),'LineWidth',2);
hold on;
plot(out_date,cfsr_wind_speed(1,:),'LineWidth',1.5);
hold on;
plot(out_date,narr_wind_speed(:,1),'LineWidth',1.5);
ylim([0 15]);xlim(datenum(2022,8:9,1));
set(gcf,'color','w');
legend({'Buoy','CFSv2','NARR'},'location','northwest');
title('Station CHBV2 - 8638901');
set(gca,'xtick',datenum(2022,8,1:5:31),'xticklabel',datestr(datenum(2022,8,1:5:31),'mm/dd'));

nexttile;
plot(out_date,buoy_wind(:,2),'LineWidth',3);
hold on;
plot(out_date,cfsr_wind_speed(2,:),'LineWidth',1.5);
hold on;
plot(out_date,narr_wind_speed(:,2),'LineWidth',1.5);
ylim([0 15]);xlim(datenum(2022,8:9,1));
set(gcf,'color','w');
title('Station 44014');
set(gca,'xtick',datenum(2022,8,1:5:31),'xticklabel',datestr(datenum(2022,8,1:5:31),'mm/dd'));
exportgraphics(gcf,'Buoy_CFSv2_NARR_WSpd.png','Resolution',300);

figure('Units','pixels','Position',[100 100 800 500])
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');


nexttile;
plot(out_date,buoy_dir(:,1),'LineWidth',2);
hold on;
plot(out_date,cfsr_dir(1,:),'LineWidth',1.5);
hold on;
plot(out_date,narr_dir(:,1),'LineWidth',1.5);
ylim([0 360]);xlim(datenum(2022,8:9,1));
set(gcf,'color','w');
title('Station CHBV2 - 8638901');
set(gca,'xtick',datenum(2022,8,1:5:31),'xticklabel',datestr(datenum(2022,8,1:5:31),'mm/dd'));

nexttile;
plot(out_date,buoy_dir(:,2),'LineWidth',3);
hold on;
plot(out_date,cfsr_dir(2,:),'LineWidth',1.5);
hold on;
plot(out_date,narr_dir(:,2),'LineWidth',1.5);
ylim([0 360]);xlim(datenum(2022,8:9,1));
set(gcf,'color','w');
legend({'Buoy','CFSv2','NARR'},'location','southwest');
title('Station 44014');
set(gca,'xtick',datenum(2022,8,1:5:31),'xticklabel',datestr(datenum(2022,8,1:5:31),'mm/dd'));
exportgraphics(gcf,'Buoy_CFSv2_NARR_WDir.png','Resolution',300);