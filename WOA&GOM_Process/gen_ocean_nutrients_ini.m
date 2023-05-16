% By Yuren Chen, 2021-06-30, Version 1.0
%--------------------------------------------------------------------------
% Program generate non-tidal boudaries from  output files of 
% GOM/WOA18 DATASET.
% Due to the mismatch of land masks between hycom
% dataset and ROMS model grid, extrapolation are performed here. 
%--------------------------------------------------------------------------

clear all; close all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
delete('./*ini.mat');

% toolbox path
addpath(path,'../../COAWST/Tools/mfiles/rutgers/utility');
addpath(path,'../../COAWST/Tools/mfiles/roms_clm');

% Output time
date_out = datenum(2001,1,1,0,0,0);

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
%--------------------- NO CHANGE BELOW-------------------------------------
%--------------------- Read from model grid--------------------------------
lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
mask = ncread(fn,'mask_rho');
angle = ncread(fn,'angle');
h = ncread(fn,'h');
sc = squeeze(set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,ones(size(h)),0,0));
sc = sc.*-1;

[r,c] = size(mask);
for i = 1:r
    for j = 1:c
        sc(i,j,:) = flipud(squeeze(sc(i,j,:)));
    end
end
%[j,i] = meshgrid(1:c,1:r);
%[j2,i2] = meshgrid(0:c+1,0:r+1);
%lon2 = griddata(i,j,lon,i2,j2,'nearest');
%lat2 = griddata(i,j,lat,i2,j2,'nearest');

%salinity
unit_c = 1;
fdir = './salinity/';
[ dat_out ] = read_woa_ini( fdir, lon, lat, mask, sc, h, date_out,'s_an');
sal = dat_out.*unit_c;

clear dat_out

%do
unit_c = 32/1025;
fdir = './do/';
[ dat_out ] = read_woa_ini( fdir, lon, lat, mask, sc, h, date_out,'o_an');
do = dat_out.*unit_c;


clear dat_out

%no3
unit_c = 14/1025;
fdir = './nitrate/';
[ dat_out ] = read_woa_ini( fdir, lon, lat, mask, sc, h, date_out,'n_an');
no3 = dat_out.*unit_c;

clear dat_out

%po4
unit_c = 31/1025;
fdir = './phosphate/';

% %!!!!!!!!!!!!!!!!!!!!!!!
% s_el_final(:,:) = 0;
% w_el_final(:,:) = 0;
% e_el_final(:,:) = 0;
% s_u_final(:,:) = 0;
% w_u_final(:,:) = 0;
% e_u_final(:,:) = 0;
% s_v_final(:,:) = 0;
% w_v_final(:,:) = 0;
% e_v_final(:,:) = 0;
% %!!!!!!!!!!!!!!!!!!!!!!!
[ dat_out ] = read_woa_ini( fdir, lon, lat, mask, sc, h, date_out,'p_an');
po4 = dat_out.*unit_c;

clear dat_out

%si
unit_c = 28/1025;
fdir = './silicate/';
[ dat_out ] = read_woa_ini( fdir, lon, lat, mask, sc, h, date_out,'i_an');
si = dat_out.*unit_c;

clear dat_out

save('ocean_nutrients_ini.mat','date_out');
save('ocean_nutrients_ini.mat','do','no3','po4','sal','si','-append','-v7.3');



