clear all;close all;

fn1 = 'WFS_2001_bry.nc';
fn2= 'WFS_2001_clm1.nc';

% toolbox path
addpath(path,'../../COAWST/Tools/mfiles/rutgers/utility');
addpath(path,'../../COAWST/Tools/mfiles/roms_clm');

% Output time
date_clm_start = datenum(2002,1,1,0,0,0);

% model grid
fn = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';

% params.
n_hycomlayer = 40;
year = 2002;
N= length(ncread(fn,'Cs_r'));
n_layer = N;
Vtransform = ncread(fn,'Vtransform');                         
Vstretching = ncread(fn,'Vstretching');
THETA_S = ncread(fn,'theta_s');                     
THETA_B = ncread(fn,'theta_b');                     
TCLINE = ncread(fn,'Tcline');
hc = ncread(fn,'hc');
interp_scheme = 1; %1 for nearest, more in the future

%list for included hycom files
hycom_dir = '../HYCOM/2002/';
fname = dir([hycom_dir 'HYCOM_*.nc']);
for i = 1:size(fname,1)
    filename{i,1} = strcat(hycom_dir,fname(i).name);
end
[n_days,c] = size(filename);
%--------------------- NO CHANGE BELOW-------------------------------------
%--------------------- Read from model grid--------------------------------
lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
lon_u = ncread(fn,'lon_u');
lon_v = ncread(fn,'lon_v');
lat_u = ncread(fn,'lat_u');
lat_v = ncread(fn,'lat_v');
mask = ncread(fn,'mask_rho');
angle = ncread(fn,'angle');
h = ncread(fn,'h');
angle_u = rho2u_2d_mw(angle);
h_u = rho2u_2d_mw(h);

[r,c] = size(mask);

% Calculate Sigma for each layer
sc = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,ones(size(h)),0,0);
sw = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,5,ones(size(h)),0,0);
sc = sc.*-1;
sw = sw.*-1;
dw = sw(:,:,1:end-1)-sw(:,:,2:end);

day = 1;
depth = double(ncread(filename{day},'depth'));
lat_hycom0 = double(ncread(filename{day},'lat'));
lon_hycom0 = double(ncread(filename{day},'lon'));

[lat_hycom,lon_hycom] = meshgrid(lat_hycom0,lon_hycom0);
[r_hycom,c_hycom] = size(lat_hycom);

x0 = repmat(lon_hycom,1,1,n_hycomlayer);
y0 = repmat(lat_hycom,1,1,n_hycomlayer);
z0 = repmat(reshape(depth,1,1,n_hycomlayer),r_hycom,c_hycom,1);

x = lon_u(200,1);
y = lat_u(200,1);
z_all = -1*set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,h(200,1),0,0);
z = z_all(end);

%--------------------- Read from model grid--------------------------------

%-------------------------3d interplation--------------------------
%
%
%temp0=zeros(r,c,n_layer,n_days);
%s0=zeros(r,c,n_layer,n_days);
%u0=zeros(r,c,n_layer,n_days);
%v0=zeros(r,c,n_layer,n_days);
date_clm = [];
k=0;
for day=1:5
    day
    tmp = double(ncread(filename{day},'time')/24+datenum(2000,01,01))-date_clm_start;
    date_clm = [date_clm;tmp];
    date_clm_m(day) = mean(double(ncread(filename{day},'time'))/24+datenum(2000,01,01))-date_clm_start;
    %temp = double(ncread(filename{day},'temperature'));
    u = double(ncread(filename{day},'u'));
    v = double(ncread(filename{day},'v'));
    u_m = mean(u,4);
    v_m = mean(v,4);
    %s = double(ncread(filename{day},'salinity'));
    for ii=1:size(u,4)
        ii
        k=k+1;
        %temp_t = temp(:,:,:,ii);
        %s_t = s(:,:,:,ii);
        u_t = u(:,:,:,ii);
        v_t = v(:,:,:,ii);

        %tmp = griddata(x0(~isnan(temp_t)),y0(~isnan(temp_t)),z0(~isnan(temp_t))...
        %    ,temp_t(~isnan(temp_t)),x,y,z);
        %tmp(isnan(tmp)) = griddata(x0(~isnan(temp_t)),y0(~isnan(temp_t)),z0(~isnan(temp_t))...
        %    ,temp_t(~isnan(temp_t)),x(isnan(tmp)),y(isnan(tmp)),z(isnan(tmp)),'nearest');
        %temp0(k) = tmp;



        %tmp = griddata(x0(~isnan(s_t)),y0(~isnan(s_t)),z0(~isnan(s_t)),...
        %    s_t(~isnan(s_t)),x,y,z,'linear');
        %tmp(isnan(tmp)) = griddata(x0(~isnan(s_t)),y0(~isnan(s_t)),z0(~isnan(s_t)),...
        %    s_t(~isnan(s_t)),x(isnan(tmp)),y(isnan(tmp)),z(isnan(tmp)),'nearest');
        %s0(k) = tmp;

        u0(k) = griddata(x0(~isnan(u_t)),y0(~isnan(u_t)),z0(~isnan(u_t)),u_t(~isnan(u_t)),x,y,z,'linear');
        v0(k) = griddata(x0(~isnan(v_t)),y0(~isnan(v_t)),z0(~isnan(v_t)),v_t(~isnan(v_t)),x,y,z,'linear');

        u_hycom(k) = cos(angle_u(200,1)).*u0(k)+sin(angle_u(200,1)).*v0(k);
        v_hycom(k) = -sin(angle_u(200,1)).*u0(k)+cos(angle_u(200,1)).*v0(k);
    end
    u0_m(day) = griddata(x0(~isnan(u_m)),y0(~isnan(u_m)),z0(~isnan(u_m)),u_m(~isnan(u_m)),x,y,z,'linear');
    v0_m(day) = griddata(x0(~isnan(v_m)),y0(~isnan(v_m)),z0(~isnan(v_m)),v_m(~isnan(v_m)),x,y,z,'linear');

    u_hycom_m(day) = cos(angle_u(200,1)).*u0_m(day)+sin(angle_u(200,1)).*v0_m(day);
    v_hycom_m(day) = -sin(angle_u(200,1)).*u0_m(day)+cos(angle_u(200,1)).*v0_m(day);
end

t1 = ncread(fn1,'v3d_time');
v1 = squeeze(ncread(fn1,'v_south',[200 1 1],[1 Inf Inf]));

t2 = ncread(fn2,'v3d_time');
v2 = squeeze(ncread(fn2,'v',[200 1 1 1],[1 1 Inf Inf]));

load('C:\Users\cheny\Desktop\EcoHAB\Non_tidal_component_preprocessing\HYCOM\non_tidal_bnd_2002.mat')

figure;
plot(date_clm,v_hycom,date_clm_m,v_hycom_m,'LineWidth',2);
%plot(date_clm,v_hycom,'LineWidth',2);
hold on;
plot(t1,v1(end,:),t2,v2(end,:));
hold on;
% for i=1:length(date_clm)/8
%     scatter(mean(date_clm([1:8]+(i-1)*8)),mean(v_hycom([1:8]+(i-1)*8)),'green','filled');
% hold on;
% end
hold on;
plot(date_out2d-date_out2d(1),s_v2d_out(200,:),'LineWidth',2);
legend('HYCOM','HYCOM_MEAN','BND','CLM','2D');
xlim([0 10]);

figure;
plot(t1,v1(end,:),t2,v2(end,:));
hold on;
plot(t1,v1(1,:),t2,v2(1,:));
hold on;
hold on;
plot(date_out2d-date_out2d(1),s_v2d_out(200,:),'LineWidth',2);
legend('BND-SUR','CLM-SUR','BND-BOT','CLM-BOT','2D');
xlim([0 10]);

figure;
plot(t1,mean(v1,1),t2,mean(v2,1),'LineWidth',4);
hold on;
plot(date_out2d-date_out2d(1),s_v2d_out(200,:),'LineWidth',2);
legend('BND','CLM','2D');
xlim([0 10]);