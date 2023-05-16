clear all; close all;

%Toolbox
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');
addpath(path,'../River_preprocessing/USGS_NWIS');

%Params.
grd = '../Model_grid/ROMS_WFS_Piney.nc';
bay_mask = ncread('../Model_grid/ROMS_WFS_Piney_bay_mask.nc','mask_rho');
%grd = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
%bay_mask = ncread('../Model_grid/ROMS_WFS_10river_grid_bay_mask.nc','mask_rho');
out_date = datenum(2021,1,1,0,0,0):1/24:datenum(2021,12,31,24,0,0); %
year = datevec(out_date(1));
year = year(1);
n_river = 11;
rname = {'STEINHATCHEE','SUWANNEE','WITHLACOOCHEE','HILLSBOROUGH',...
    'ALAFIA','LITTLE_MANATEE','MANATEE','MYAKKA','PEACE','CALOOSAHATCHEE','LAKE'};
varname = {'NO23','NH4','ON','TN','PO4','TP','SI','OC','OP'};
on_flag = [0 0 1 0 0 0 1 1 1 1 1];  %River ON calculation flag
include_cruise = [0 0 0 0];%NO23, PO4, SI, DON 
dep_range = [0 1000];  %depth range for selected cruise data
time_range = [datenum(year-1,9,1) datenum(year,4,1)];%time range for selected IN cruise data
time_range2 = [datenum(year,1,1) datenum(year,5,1)];%time range for selected DON cruise data

skip_step = [1 1 1 1 1 0];

% Flag for using LOADEST
%           'NO23' 'NH4' 'ON' 'TN' 'PO4' 'TP' 'SI' 'OC' 'OP'             
var_flag(1,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %STEINHATCHEE
var_flag(2,:) = ...
          [    0     0    0    0     0    0    0    0    0]; %SUWANNEE
var_flag(3,:) = ...
          [    1     1    1    1     1    1    0    1    1]; %WITHLACOOCHEE
var_flag(4,:) = ...
          [    0     0    0    0     0    0    1    0    0]; %HILLSBOROUGH
var_flag(5,:) = ...
          [    0     0    0    0     0    0    1    0    0]; %ALAFIA
var_flag(6,:) = ...
          [    0     0    0    0     0    0    1    0    0]; %LITTLE_MANATEE
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
[nutri_dat,nutri_mean] = pre_process_river_nutrients('2021-high.xlsx',varname,n_river,on_flag); %!!!!Unfinished
[temp_dat,temp_mean] = pre_preocess_river_temp('River_temp.xlsx',n_river);
[sal_dat,sal_mean] = pre_preocess_river_salt('River_sal.xlsx',n_river);

save(strcat('WA_',num2str(year),'.mat'),'nutri_mean','nutri_dat','temp_mean','temp_dat','sal_mean','sal_dat');
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

save(strcat('calib_',num2str(year),'.mat'),'calflow','caldat');
end

%Step 3 run loadest 
if(skip_step(3)==0)
disp('Please run loadest manually! Press any key after finishing...');
pause;
end

%Step 4 read_loadest
if(skip_step(4)==0)
[result_flow,result_conc] = ...
    read_plot_loadest_result(rname, varname, var_flag, n_river, out_date, './fig/');

save(strcat('loadest_',num2str(year),'.mat'),'result_conc','result_flow');
end
%Step 5 river boundaries
if(skip_step(5)==0)
gen_river_nutrients(grd, varname, var_flag, n_river, out_date);
end

%Step 6 bay initial conditions
if(skip_step(6)==0)
gen_bay_nutrients_ini('../Cruise_data',['./bay_initials_',num2str(year),'.xlsx']...
    ,grd,bay_mask,year,include_cruise,dep_range,time_range,time_range2);
end