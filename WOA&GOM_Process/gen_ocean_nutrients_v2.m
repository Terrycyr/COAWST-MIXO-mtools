 % By Yuren Chen, 2021-11-7, Version 2.0
%--------------------------------------------------------------------------
% Program generate non-tidal boudaries from  output files of 
% GOM/WOA18 DATASET.
% Due to the mismatch of land masks between hycom
% dataset and ROMS model grid, extrapolation are performed here. 
%--------------------------------------------------------------------------

clear all; close all;

addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

% toolbox path
addpath(path,'../../COAWST/Tools/mfiles/rutgers/utility');
addpath(path,'../../COAWST/Tools/mfiles/roms_clm');

% Output time
date_out = datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,24,0,0);

% model grid
fn = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';

% params.
year = 2001;
N= length(ncread(fn,'Cs_r'));
Vtransform = ncread(fn,'Vtransform');                         
Vstretching = ncread(fn,'Vstretching');
THETA_S = ncread(fn,'theta_s');                     
THETA_B = ncread(fn,'theta_b');                     
TCLINE = ncread(fn,'Tcline');
hc = ncread(fn,'hc');
dir_flag = [1 1 0 1]; %S N W E !!!!!!!!!!!!!!!!!!!!
%--------------------- NO CHANGE BELOW-------------------------------------
%--------------------- Read from model grid--------------------------------
lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
mask = ncread(fn,'mask_rho');
angle = ncread(fn,'angle');
h = ncread(fn,'h');
[r,c] = size(mask);

load(strcat('ocean_nutrients_ini_',num2str(year),'.mat'));
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
no3 = no3_raw;
po4 = po4_raw;
si = si_raw;

s_lon = lon(:,1);
s_lat = lat(:,1);
s_h = h(:,1);
s_mask = mask(:,1);
s_angle = angle(:,1);

n_lon = lon(:,end);
n_lat = lat(:,end);
n_h = h(:,end);
n_mask = mask(:,end);
n_angle = angle(:,end);

w_lon = lon(1,:);
w_lat = lat(1,:);
w_h = h(1,:);
w_mask = mask(1,:);
w_angle = angle(1,:);


e_lon = lon(end,:);
e_lat = lat(end,:);
e_h = h(end,:);
e_mask = mask(end,:);
e_angle = angle(end,:);

sb_num = length(s_lon(:,1));
nb_num = length(n_lon(:,1));
wb_num = length(w_lon(1,:));
eb_num = length(e_lon(1,:));

sb_range = [1:sb_num];
nb_range = [1:nb_num]+sb_range(end);
eb_range = [1:eb_num]+nb_range(end);
wb_range = [1:wb_num]+eb_range(end);


% Calculate Sigma for each layer
sc = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,ones(size(h)),0,0);
sw = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,5,ones(size(h)),0,0);
sc = sc.*-1;
sw = sw.*-1;
dw = sw(:,:,1:end-1)-sw(:,:,2:end);
s_sc = squeeze(sc(:,1,:));
n_sc = squeeze(sc(:,end,:));
e_sc = squeeze(sc(end,:,:));
w_sc = squeeze(sc(1,:,:));

s_dw = squeeze(dw(:,1,:));
n_dw = squeeze(dw(:,end,:));
e_dw = squeeze(dw(end,:,:));
w_dw = squeeze(dw(1,:,:));

% all boundary points
b_lon = [s_lon(:,1)',n_lon(:,1)',e_lon(1,:),w_lon(1,:)];
b_lat = [s_lat(:,1)',n_lat(:,1)',e_lat(1,:),w_lat(1,:)];
b_sc = [s_sc;n_sc;e_sc;w_sc];
b_sc = flipud(b_sc')';
b_h = [s_h(:,1)',n_h(:,1)',e_h(1,:),w_h(1,:)];
b_mask = [s_mask(:,1)',n_mask(:,1)',e_mask(1,:),w_mask(1,:)];
b_angle = [s_angle(:,1)',n_angle(:,1)',e_angle(1,:),w_angle(1,:)];
b_num = length(b_lon);


%salinity
w_sal = squeeze(sal(1,:,:));
e_sal = squeeze(sal(end,:,:));
n_sal = squeeze(sal(:,end,:));
s_sal = squeeze(sal(:,1,:));

w_sal = repmat(w_sal,1,1,length(date_out));
e_sal = repmat(e_sal,1,1,length(date_out));
n_sal = repmat(n_sal,1,1,length(date_out));
s_sal = repmat(s_sal,1,1,length(date_out));

%do
w_do = squeeze(do(1,:,:));
e_do = squeeze(do(end,:,:));
n_do = squeeze(do(:,end,:));
s_do = squeeze(do(:,1,:));
w_do = repmat(w_do,1,1,length(date_out));
e_do = repmat(e_do,1,1,length(date_out));
n_do = repmat(n_do,1,1,length(date_out));
s_do = repmat(s_do,1,1,length(date_out));

%no3
w_no3 = squeeze(no3(1,:,:));
e_no3 = squeeze(no3(end,:,:));
n_no3 = squeeze(no3(:,end,:));
s_no3 = squeeze(no3(:,1,:));

w_no3 = repmat(w_no3,1,1,length(date_out));
e_no3 = repmat(e_no3,1,1,length(date_out));
n_no3 = repmat(n_no3,1,1,length(date_out));
s_no3 = repmat(s_no3,1,1,length(date_out));

%nh4
w_nh4 = squeeze(nh4(1,:,:));
e_nh4 = squeeze(nh4(end,:,:));
n_nh4 = squeeze(nh4(:,end,:));
s_nh4 = squeeze(nh4(:,1,:));

w_nh4 = repmat(w_nh4,1,1,length(date_out));
e_nh4 = repmat(e_nh4,1,1,length(date_out));
n_nh4 = repmat(n_nh4,1,1,length(date_out));
s_nh4 = repmat(s_nh4,1,1,length(date_out));

%po4
w_po4 = squeeze(po4(1,:,:));
e_po4 = squeeze(po4(end,:,:));
n_po4 = squeeze(po4(:,end,:));
s_po4 = squeeze(po4(:,1,:));

w_po4 = repmat(w_po4,1,1,length(date_out));
e_po4 = repmat(e_po4,1,1,length(date_out));
n_po4 = repmat(n_po4,1,1,length(date_out));
s_po4 = repmat(s_po4,1,1,length(date_out));

%si
w_si = squeeze(si(1,:,:));
e_si = squeeze(si(end,:,:));
n_si = squeeze(si(:,end,:));
s_si = squeeze(si(:,1,:));

w_si = repmat(w_si,1,1,length(date_out));
e_si = repmat(e_si,1,1,length(date_out));
n_si = repmat(n_si,1,1,length(date_out));
s_si = repmat(s_si,1,1,length(date_out));

%temp
w_temp = squeeze(temp(1,:,:));
e_temp = squeeze(temp(end,:,:));
n_temp = squeeze(temp(:,end,:));
s_temp = squeeze(temp(:,1,:));

w_temp = repmat(w_temp,1,1,length(date_out));
e_temp = repmat(e_temp,1,1,length(date_out));
n_temp = repmat(n_temp,1,1,length(date_out));
s_temp = repmat(s_temp,1,1,length(date_out));

clear dat_out

save(strcat('ocean_nutrients_bnd_',num2str(year),'.mat'),'date_out');
if(dir_flag(1)==1)
    save(strcat('ocean_nutrients_bnd_',num2str(year),'.mat'),'s_do','s_no3','s_nh4','s_po4','s_sal','s_si','s_temp','-append','-v7.3');
end
if(dir_flag(2)==1)
    save(strcat('ocean_nutrients_bnd_',num2str(year),'.mat'),'n_do','n_no3','n_nh4','n_po4','n_sal','n_si','n_temp','-append','-v7.3');
end
if(dir_flag(3)==1)
    save(strcat('ocean_nutrients_bnd_',num2str(year),'.mat'),'w_do','w_no3','w_nh4','w_po4','w_sal','w_si','w_temp','-append','-v7.3');
end
if(dir_flag(4)==1)
    save(strcat('ocean_nutrients_bnd_',num2str(year),'.mat'),'e_do','e_no3','e_nh4','e_po4','e_sal','e_si','e_temp','-append','-v7.3');
end

