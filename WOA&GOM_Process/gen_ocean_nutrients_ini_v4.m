% By Yuren Chen, 2021-11-19, Version 3.1
%--------------------------------------------------------------------------
% Program generate initial fields from  output files of 
% GOM/WOA18/NEGOM DATASET.
% For nutrients, mean profiles are used instead of spatial interpolation. 
%--------------------------------------------------------------------------

clear all; close all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

% toolbox path
addpath(path,'../../COAWST/Tools/mfiles/rutgers/utility');
addpath(path,'../../COAWST/Tools/mfiles/roms_clm');

% Output time
date_out = datenum(2021,3,30,0,0,0);

% model grid
%fn = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
fn = '../Model_grid/ROMS_WFS_Piney.nc';

% params.
year = 2021;
N= length(ncread(fn,'Cs_r'));
Vtransform = ncread(fn,'Vtransform');                         
Vstretching = ncread(fn,'Vstretching');
THETA_S = ncread(fn,'theta_s');                     
THETA_B = ncread(fn,'theta_b');                     
TCLINE = ncread(fn,'Tcline');
hc = ncread(fn,'hc');
%--------------------- NO CHANGE BELOW-------------------------------------
%--------------------- Read from model grid--------------------------------
delete(strcat('ocean_nutrients_ini_',num2str(year),'.mat'));
lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
mask = ncread(fn,'mask_rho');
angle = ncread(fn,'angle');
h = ncread(fn,'h');
sc = squeeze(set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,ones(size(h)),0,0));
sc = sc.*-1;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tmp = h(412:end,1:end);
tmp(tmp>60) = 0;
h2 = h;
h2(412:end,1:end) = tmp;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
fdir = '../GOM18/salinity/';
[ dat_out ] = read_woa_ini( fdir, lon, lat, sc, h, date_out,'s_an');
sal = dat_out.*unit_c;

clear dat_out

%do
unit_c = 32/1025;
fdir = '../GOM18/do/';
[ dat_out ] = read_woa_ini3( fdir, lon, lat, sc, h, date_out,'o_an',0);
do = dat_out.*unit_c;

clear dat_out

%no3
unit_c = 14/1025;
%fdir = '../GOM18/nitrate/';
[ dat_out ] = read_negom_ini( '../GOM18/NEGOM_nutrients.xlsx', lon, lat, sc, h, 4);
no3 = dat_out.*unit_c;
clear dat_out

%nh4
unit_c = 14/1025;
[ dat_out ] = read_negom_ini( '../GOM18/NEGOM_nutrients.xlsx', lon, lat, sc, h, 3);
nh4 = dat_out.*unit_c;

clear dat_out

%DON
unit_c = 14/1025;
[ dat_out ] = read_ecohab_ini( '../Cruise_data/DON_BOT.xlsx', lon, lat, sc, h, 1);
don = dat_out.*unit_c;

clear dat_out

%po4
unit_c = 31/1025;
fdir = '../GOM18/phosphate/';
[ dat_out ] = read_negom_ini( '../GOM18/NEGOM_nutrients.xlsx', lon, lat, sc, h, 1);
po4 = dat_out.*unit_c;

clear dat_out

%si
unit_c = 28/1025;
fdir = '../GOM18/silicate/';
[ dat_out, si_m, si_dep] = read_woa_ini3( fdir, lon, lat, sc, h, date_out,'i_an',0);
si = dat_out.*unit_c;
clear dat_out

%temp
unit_c = 1;
fdir = '../GOM18/temp/';
[ dat_out ] = read_woa_ini( fdir, lon, lat, sc, h, date_out,'t_an');
temp = dat_out.*unit_c;

clear dat_out

%for ROMS
[r,c] = size(mask);
for i = 1:r
    for j = 1:c
        no3_roms(i,j,:) = flipud(squeeze(no3(i,j,:)));
        sal_roms(i,j,:) = flipud(squeeze(sal(i,j,:)));
        temp_roms(i,j,:) = flipud(squeeze(temp(i,j,:)));
    end
end

date_ini = date_out;

save(strcat('ocean_nutrients_ini_',num2str(year),'.mat'),'date_ini');
save(strcat('ocean_nutrients_ini_',num2str(year),'.mat'),'do','no3','nh4','po4','sal','si','temp','don','-append','-v7.3');
%save(strcat('ocean_nutrients_ini_',num2str(year),'.mat'),'no3_raw','po4_raw','si_raw','-append','-v7.3');
save(strcat('ocean_nutrients_ini_',num2str(year),'.mat'),'no3_roms','sal_roms','temp_roms','-append','-v7.3');
