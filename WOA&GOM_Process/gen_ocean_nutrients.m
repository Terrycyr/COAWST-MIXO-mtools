% By Yuren Chen, 2021-06-30, Version 1.0
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
fn = '../Model_grid/ROMS_WFS_10river_grid_v10.nc';

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
unit_c = 1;
fdir = './salinity/';
[ dat_out ] = read_woa( fdir, b_lon, b_lat, b_mask, b_sc, b_h, date_out,'s_an');
s_sal = dat_out(sb_range,:,:).*unit_c;
n_sal = dat_out(nb_range,:,:).*unit_c;
e_sal = dat_out(eb_range,:,:).*unit_c;
w_sal = dat_out(wb_range,:,:).*unit_c;

s_sal = extend_mat(s_sal,1);
n_sal = extend_mat(n_sal,1);
e_sal = extend_mat(e_sal,1);
w_sal = extend_mat(w_sal,1);

clear dat_out

%do
unit_c = 32/1025;
fdir = './do/';
[ dat_out ] = read_woa2( fdir, b_lon, b_lat, b_mask, b_sc, b_h, date_out,'o_an');
s_do = dat_out(sb_range,:,:).*unit_c;
n_do = dat_out(nb_range,:,:).*unit_c;
e_do = dat_out(eb_range,:,:).*unit_c;
w_do = dat_out(wb_range,:,:).*unit_c;

s_do = extend_mat(s_do,1);
n_do = extend_mat(n_do,1);
e_do = extend_mat(e_do,1);
w_do = extend_mat(w_do,1);

clear dat_out

%no3
unit_c = 14/1025;
fdir = './nitrate/';
[ dat_out, dat_m ] = read_woa2( fdir, b_lon, b_lat, b_mask, b_sc, b_h, date_out,'n_an');
s_no3 = dat_out(sb_range,:,:).*unit_c;
n_no3 = dat_out(nb_range,:,:).*unit_c;
e_no3 = dat_out(eb_range,:,:).*unit_c;
w_no3 = dat_out(wb_range,:,:).*unit_c;

s_no3 = extend_mat(s_no3,1);
n_no3 = extend_mat(n_no3,1);
e_no3 = extend_mat(e_no3,1);
w_no3 = extend_mat(w_no3,1);

clear dat_out

%po4
unit_c = 31/1025;
fdir = './phosphate/';
[ dat_out ] = read_woa2( fdir, b_lon, b_lat, b_mask, b_sc, b_h, date_out,'p_an');
s_po4 = dat_out(sb_range,:,:).*unit_c;
n_po4 = dat_out(nb_range,:,:).*unit_c;
e_po4 = dat_out(eb_range,:,:).*unit_c;
w_po4 = dat_out(wb_range,:,:).*unit_c;

s_po4 = extend_mat(s_po4,1);
n_po4 = extend_mat(n_po4,1);
e_po4 = extend_mat(e_po4,1);
w_po4 = extend_mat(w_po4,1);

clear dat_out

%si
unit_c = 28/1025;
fdir = './silicate/';
[ dat_out ] = read_woa2( fdir, b_lon, b_lat, b_mask, b_sc, b_h, date_out,'i_an');
s_si = dat_out(sb_range,:,:).*unit_c;
n_si = dat_out(nb_range,:,:).*unit_c;
e_si = dat_out(eb_range,:,:).*unit_c;
w_si = dat_out(wb_range,:,:).*unit_c;

s_si = extend_mat(s_si,1);
n_si = extend_mat(n_si,1);
e_si = extend_mat(e_si,1);
w_si = extend_mat(w_si,1);

%temp
unit_c = 1;
fdir = './temp/';
[ dat_out ] = read_woa( fdir, b_lon, b_lat, b_mask, b_sc, b_h, date_out,'t_an');
s_temp = dat_out(sb_range,:,:).*unit_c;
n_temp = dat_out(nb_range,:,:).*unit_c;
e_temp = dat_out(eb_range,:,:).*unit_c;
w_temp = dat_out(wb_range,:,:).*unit_c;

s_temp = extend_mat(s_temp,1);
n_temp = extend_mat(n_temp,1);
e_temp = extend_mat(e_temp,1);
w_temp = extend_mat(w_temp,1);

clear dat_out

save('ocean_nutrients_bnd.mat','date_out');
if(dir_flag(1)==1)
    save('ocean_nutrients_bnd.mat','s_do','s_no3','s_po4','s_sal','s_si','s_temp','-append','-v7.3');
end
if(dir_flag(2)==1)
    save('ocean_nutrients_bnd.mat','n_do','n_no3','n_po4','n_sal','n_si','n_temp','-append','-v7.3');
end
if(dir_flag(3)==1)
    save('ocean_nutrients_bnd.mat','w_do','w_no3','w_po4','w_sal','w_si','w_temp','-append','-v7.3');
end
if(dir_flag(4)==1)
    save('ocean_nutrients_bnd.mat','e_do','e_no3','e_po4','e_sal','e_si','e_temp','-append','-v7.3');
end



