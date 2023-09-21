% By Yuren Chen, 2021-04-06, Version 2.0
% Fixed bugs, 2021-06-18, Version 2.1
% Change to 3d interpolation, 2022-01-24. Version 2.2
%--------------------------------------------------------------------------
% Program generate non-tidal boudaries from  output files of 
% Simple Ocean Data Assimilation (hycom).
% Due to the mismatch of land masks between hycom
% dataset and ROMS model grid, extrapolation are performed here. 
%--------------------------------------------------------------------------

clear all; close all;

% toolbox path
addpath(path,'../../../COAWST/Tools/mfiles/rutgers/utility');
addpath(path,'../../../COAWST/Tools/mfiles/roms_clm');

% Output time
year = 2022;
date_out3d = datenum(year,1,1,0,0,0):6/24:datenum(year,12,31,24,0,0);
date_out2d = datenum(year,1,1,0,0,0):1/24:datenum(year,12,31,24,0,0);

% model grid  !!!!!!!!!!!!!
%fn = '../../Model_grid/ROMS_WFS_10river_grid_v11.nc';
fn = '../../Model_grid/ROMS_WFS_piney.nc';

% params.
n_hycomlayer = 40;
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
hycom_dir = ['../../HYCOM/',num2str(year),'/'];
fname = dir([hycom_dir 'HYCOM_*.nc']);
filename{1,1} = strcat(hycom_dir,fname(1).name);

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

[r,c] = size(mask);

% Calculate Sigma for each layer
sc = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,ones(size(h)),0,0);
sw = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,5,ones(size(h)),0,0);
sc = sc.*-1;
sw = sw.*-1;
dw = sw(:,:,1:end-1)-sw(:,:,2:end);

day = 1;
depth = double(ncread(filename{day},'depth'));
lat_hycom = double(ncread(filename{day},'lat'));
lon_hycom = double(ncread(filename{day},'lon'));
if(sum(size(lat_hycom)>1)==1)
    [lat_hycom,lon_hycom] = meshgrid(lat_hycom,lon_hycom);
end
[r_hycom,c_hycom] = size(lat_hycom);

x0 = repmat(lon_hycom,1,1,n_hycomlayer);
y0 = repmat(lat_hycom,1,1,n_hycomlayer);
z0 = repmat(reshape(depth,1,1,n_hycomlayer),r_hycom,c_hycom,1);
x = repmat(lon,1,1,n_layer);
y = repmat(lat,1,1,n_layer);
z = -1*set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,h,0,0);

%--------------------- Read from model grid--------------------------------

%-------------------------horizontal interplation--------------------------
% Barnes method is used in the horizontal interpolation and extrapolation
% of temp., sal. and zeta, while that of u and v is linear interoplation to
% avoid introducing extra flux. 
data_time = [];
for day=1:1
    day
    tmp = double(ncread(filename{day},'time'))/24+datenum(2000,01,01);
    data_time = [data_time;tmp];   
    temp = double(ncread(filename{day},'temperature'));
    u = double(ncread(filename{day},'u'));
    v = double(ncread(filename{day},'v'));
    s = double(ncread(filename{day},'salinity'));
    el = double(ncread(filename{day},'surf_el'));
   
    el_t = el(:,:,1);
    
%    el0(:,i+nrec_total) = barnes(lon_hycom(~isnan(el_t)),lat_hycom(~isnan(el_t)),el_t(~isnan(el_t)),lon,lat,1,1,3);
    el0(:,:) = griddata(lon_hycom(~isnan(el_t)),lat_hycom(~isnan(el_t)),el_t(~isnan(el_t)),lon,lat,'linear');
    
    temp_t = mean(temp,4);
    s_t = mean(s,4);
    u_t = mean(u,4);
    v_t = mean(v,4);

    tmp = griddata(x0(~isnan(temp_t)),y0(~isnan(temp_t)),z0(~isnan(temp_t))...
        ,temp_t(~isnan(temp_t)),x,y,z);
    tmp(isnan(tmp)) = griddata(x0(~isnan(temp_t)),y0(~isnan(temp_t)),z0(~isnan(temp_t))...
        ,temp_t(~isnan(temp_t)),x(isnan(tmp)),y(isnan(tmp)),z(isnan(tmp)),'nearest');

    temp0(:,:,:,day) = tmp;



    tmp = griddata(x0(~isnan(s_t)),y0(~isnan(s_t)),z0(~isnan(s_t)),...
        s_t(~isnan(s_t)),x,y,z,'linear');
    tmp(isnan(tmp)) = griddata(x0(~isnan(s_t)),y0(~isnan(s_t)),z0(~isnan(s_t)),...
        s_t(~isnan(s_t)),x(isnan(tmp)),y(isnan(tmp)),z(isnan(tmp)),'nearest');

    s0(:,:,:,day) = tmp;

    u0(:,:,:,day) = griddata(x0(~isnan(u_t)),y0(~isnan(u_t)),z0(~isnan(u_t)),u_t(~isnan(u_t)),x,y,z,'linear');
    v0(:,:,:,day) = griddata(x0(~isnan(v_t)),y0(~isnan(v_t)),z0(~isnan(v_t)),v_t(~isnan(v_t)),x,y,z,'linear');
end

el0(isnan(el0)) = 0;
u0(isnan(u0)) = 0;
v0(isnan(v0)) = 0;


temp3=temp0;
s3=s0;
u3=u0;
v3=v0;
el3=el0;

%--------------------------rotate u,v to x,e-------------------------------
%rotate u,v to x,e

h_u = rho2u_2d_mw(h);
h_v = rho2v_2d_mw(h);
angle_u = rho2u_2d_mw(angle);
angle_v = rho2v_2d_mw(angle);
mask_u = rho2u_2d_mw(mask);
mask_v = rho2v_2d_mw(mask);

for n_l= 1:n_layer
    u4(:,:,n_l) = rho2u_2d_mw(squeeze(u3(:,:,n_l)));
    v4(:,:,n_l) = rho2v_2d_mw(squeeze(v3(:,:,n_l)));
    uv3(:,:,n_l) = rho2u_2d_mw(squeeze(v3(:,:,n_l)));
    vu3(:,:,n_l) = rho2v_2d_mw(squeeze(u3(:,:,n_l)));
    
    u5(:,:,n_l) = cos(angle_u).*u4(:,:,n_l)+sin(angle_u).*uv3(:,:,n_l);
    v5(:,:,n_l) = -sin(angle_v).*vu3(:,:,n_l)+cos(angle_v).*v4(:,:,n_l);
    
    tmp = u5(:,:,n_l);
    tmp(mask_u==0) = 0;
    u5(:,:,n_l) = tmp;
    
    tmp = v5(:,:,n_l);
    tmp(mask_v==0) = 0;
    v5(:,:,n_l) = tmp;
    
    tmp = temp3(:,:,n_l);
    tmp(mask==0) = 0;
    temp3(:,:,n_l) = tmp;
    
    tmp = s3(:,:,n_l);
    tmp(mask==0) = 0;
    s3(:,:,n_l) = tmp;
end    
%--------------------------rotate u,v to x,e-------------------------------

%---------------------------------2D  VELOCITY-----------------------------
dw_u = dw(1:end-1,:,:);
dw_v = dw(:,1:end-1,:);
u2d_ini_out = sum(u5.*dw_u,3);
v2d_ini_out = sum(v5.*dw_v,3);

%Output

el_ini_out = el0;
s_ini_out = s3;
temp_ini_out = temp3;
u_ini_out = u5;
v_ini_out = v5;

save(strcat('non_tidal_ini_',num2str(year),'.mat'),'data_time','-v7.3');
save(strcat('non_tidal_ini_',num2str(year),'.mat'),'el_ini_out','-append');           
save(strcat('non_tidal_ini_',num2str(year),'.mat'),'s_ini_out','-append');
save(strcat('non_tidal_ini_',num2str(year),'.mat'),'temp_ini_out','-append');
save(strcat('non_tidal_ini_',num2str(year),'.mat'),'u_ini_out','-append');
save(strcat('non_tidal_ini_',num2str(year),'.mat'),'v_ini_out','-append');
save(strcat('non_tidal_ini_',num2str(year),'.mat'),'u2d_ini_out','-append');
save(strcat('non_tidal_ini_',num2str(year),'.mat'),'v2d_ini_out','-append');

%Figures for checking
figure;
quiver(lon_u(:,1:end-1),lat_u(:,1:end-1),u4(:,1:end-1,end),v4(1:end-1,:,end),0,'r');
hold on;
quiver(lon_u(:,1:end-1),lat_u(:,1:end-1),u_ini_out(:,1:end-1,end),v_ini_out(1:end-1,:,end),0);
title('Surface Velocity');
legend('Before rotation','After rotation');

figure;
subplot(1,2,1);
contourf(lon,lat,temp3(:,:,end));
title('Surface Temp.');
subplot(1,2,2);
contourf(lon,lat,temp3(:,:,1));
title('Bottom Temp.');


figure;
subplot(1,2,1);
contourf(lon,lat,s3(:,:,end));
title('Surface Sal.');
subplot(1,2,2);
contourf(lon,lat,s3(:,:,1));
title('Bottom Sal.');





