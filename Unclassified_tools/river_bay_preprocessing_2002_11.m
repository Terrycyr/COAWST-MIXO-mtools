clear all; close all;

%Toolbox
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');
addpath(path,'../River_preprocessing/USGS_NWIS');

%Params.
grd = '../Model_grid/ROMS_WFS_new.nc';
bay_mask = ncread('../Model_grid/ROMS_WFS_new_bay_mask.nc','mask_rho');
year = 2002;
include_cruise = [0 0 0 0];%NO23, PO4, SI, DON 
Ini_DON_method = 'TN';
bay_tss_flag = 0;

fn_bayini = strcat('bay_ini_',num2str(year),'_11.mat');

%Step 6 bay initial conditions
gen_bay_nutrients_ini('../Cruise_data',['./bay_initials_',num2str(year),'_11.xlsx']...
    ,fn_bayini,grd,bay_mask,year,include_cruise,[],[],[],Ini_DON_method,bay_tss_flag);