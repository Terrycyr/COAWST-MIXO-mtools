clear all; close all;

%Toolbox
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');
addpath(path,'../River_preprocessing/USGS_NWIS');

%Params.
grd = '../Model_grid/ROMS_WFS_new.nc';
bay_mask = ncread('../Model_grid/ROMS_WFS_new_bay_mask.nc','mask_rho');
out_date = datenum(2017,1,1,0,0,0):datenum(2017,12,31,24,0,0); %
year = datevec(out_date(1));
year = year(1);
n_river = 11;
rname = {'STEINHATCHEE','SUWANNEE','WITHLACOOCHEE','HILLSBOROUGH',...
    'ALAFIA','LITTLE_MANATEE','MANATEE','MYAKKA','PEACE','CALOOSAHATCHEE','LAKE'};
varname = {'NO23','NH4','ON','TN','PO4','TP','SI','OC','OP'};
on_flag = [1 1 1 0 0 0 1 1 1 1 1];  %River ON calculation flag
include_cruise = [0 0 0 0];%NO23, PO4, SI, DON 
Ini_DON_method = 'TN';
bay_tss_flag = 0;
skip_step = [0 0 0 0 0 1];

fn_WA = strcat('WA_',num2str(year),'.mat');
fn_calib = strcat('calib_',num2str(year),'.mat');
fn_loadest = strcat('loadest_',num2str(year),'.mat');
fn_rout = strcat('WA_river_bnd_',num2str(year),'.mat');
fn_bayini = strcat('bay_ini_',num2str(year),'.mat');
% Flag for using LOADEST
%           'NO23' 'NH4' 'ON' 'TN' 'PO4' 'TP' 'SI' 'OC' 'OP'             
var_flag(1,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %STEINHATCHEE
var_flag(2,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %SUWANNEE
var_flag(3,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %WITHLACOOCHEE
var_flag(4,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %HILLSBOROUGH
var_flag(5,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %ALAFIA
var_flag(6,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %LITTLE_MANATEE
var_flag(7,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %MANATEE
var_flag(8,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %MYAKKA
var_flag(9,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %PEACE
var_flag(10,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %CALOOSAHATCHEE
var_flag(11,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %LAKE

%Step 1
%Read river raw data
if(skip_step(1)==0)
[nutri_dat,nutri_mean] = pre_process_river_nutrients('2016-2019-high.xlsx',varname,n_river,on_flag); %!!!!Unfinished
[temp_dat,temp_mean] = pre_preocess_river_temp('River_temp.xlsx',n_river);
[sal_dat,sal_mean] = pre_preocess_river_salt('River_sal.xlsx',n_river);

save(fn_WA,'nutri_mean','nutri_dat','temp_mean','temp_dat','sal_mean','sal_dat');
end

%Step 2
%Create_loadest files
%River discharge data
if(skip_step(2)==0)
r_dir = '../USGS_River/';
r(1) = load([r_dir,'USGS_02324000_STEINHATCHEE.mat']);
r(2) = load([r_dir,'USGS_02323500_SUWANNEE.mat']);
r(3) = load([r_dir,'USGS_02313250_WITHLACOOCHEE.mat']);
r(4) = load([r_dir,'USGS_02304500_HILLSBOROUGH.mat']);
r(5) = load([r_dir,'USGS_02301500_ALAFIA.mat']);
r(6) = load([r_dir,'USGS_02300500_LITTLE_MANATEE.mat']);
r(7) = load([r_dir,'USGS_02299950_MANATEE.mat']);
r(8) = load([r_dir,'USGS_02298880_MYAKKA.mat']);
r(9) = load([r_dir,'USGS_02296750_PEACE.mat']);
r(10) = load([r_dir,'USGS_02292900_CALOOSAHATCHEE.mat']);
r(11) = load([r_dir,'USGS_02307498_LAKE_TARPON.mat']);

[calflow,caldat] = create_loadest_file(r, rname, varname, var_flag, n_river, out_date);

save(fn_calib,'calflow','caldat');
end

%Step 3 run loadest 
if(skip_step(3)==0)
disp('Please run loadest manually! Press any key after finishing...');
pause;
end

%Step 4 read_loadest
if(skip_step(4)==0)
    [result_flow,result_conc] = ...
        read_plot_loadest_result(rname, varname, var_flag, n_river, out_date...
        , './fig/',fn_WA,fn_calib);

    save(fn_loadest,'result_conc','result_flow');
end
%Step 5 river boundaries
if(skip_step(5)==0)
gen_river_nutrients(grd, varname, var_flag, n_river, out_date,fn_WA,fn_loadest,fn_rout);
end

%Step 6 bay initial conditions
if(skip_step(6)==0)
gen_bay_nutrients_ini('../Cruise_data',['./bay_initials_',num2str(year),'.xlsx']...
    ,fn_bayini,grd,bay_mask,year,include_cruise,[],[],[],Ini_DON_method,bay_tss_flag);
end