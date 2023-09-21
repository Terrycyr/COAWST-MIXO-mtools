clear all; close all;

%Toolbox
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');
addpath(path,'../River_preprocessing/USGS_NWIS');

%Params.
grd = '../Model_grid/ROMS_WFS_new.nc';
bay_mask = ncread('../Model_grid/ROMS_WFS_new_bay_mask.nc','mask_rho');
out_date = datenum(2022,1,1,0,0,0):datenum(2022,12,31,24,0,0); %
year = datevec(out_date(1));
year = year(1);
n_river = 11;
rname = {'STEINHATCHEE','SUWANNEE','WITHLACOOCHEE','HILLSBOROUGH',...
    'ALAFIA','LITTLE_MANATEE','MANATEE','MYAKKA','PEACE','CALOOSAHATCHEE','LAKE'};
varname = {'NO23','NH4','ON','TN','PO4','TP','SI','OC','OP','TSS'};
on_flag = [1 1 1 0 0 0 1 0 0 0 0];  %River ON calculation flag
include_cruise = [0 0 0 0];%NO23, PO4, SI, DON 
dep_range = [0 1000];  %depth range for selected cruise data
time_range = [datenum(year-1,9,1) datenum(year,4,1)];%time range for selected IN cruise data
time_range2 = [datenum(year,1,1) datenum(year,5,1)];%time range for selected DON cruise data
Ini_DON_method = 'TKN';
bay_tss_flag = 1;
skip_step = [0 0 0 0 0 0 0];


%Step 1
%Read river raw data
if(skip_step(1)==0)
    [nutri_dat,nutri_mean] = pre_process_river_nutrients('2022-high.xlsx',varname(1:end-1),n_river,on_flag);  %!!!!
    [temp_dat,temp_mean] = pre_preocess_river_temp('River_temp.xlsx',n_river);
    [sal_dat,sal_mean] = pre_preocess_river_salt('River_sal.xlsx',n_river);
    [tss_dat,tss_mean] = pre_preocess_river_tss('River_tss.xlsx',n_river);

    save(strcat('WA_',num2str(year),'.mat'),'nutri_mean','nutri_dat','temp_mean','temp_dat','sal_mean','sal_dat');
    save(strcat('WA_',num2str(year),'.mat'),'tss_mean','tss_dat','-append');
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
%           'NO23' 'NH4' 'ON' 'TN' 'PO4' 'TP' 'SI' 'OC' 'OP'  'TSS'        
var_flag(1,:) = ...
          [    0     0    0    0     0    0    0    0    0    1]; %STEINHATCHEE
var_flag(2,:) = ...
          [    0     0    0    0     0    0    0    0    0    1]; %SUWANNEE
var_flag(3,:) = ...
          [    1     1    1    1     1    1    0    1    1    1]; %WITHLACOOCHEE
var_flag(4,:) = ...
          [    0     0    0    0     0    0    1    0    0    1]; %HILLSBOROUGH
var_flag(5,:) = ...
          [    0     0    0    0     0    0    1    0    0    1]; %ALAFIA
var_flag(6,:) = ...
          [    0     0    0    0     0    0    1    0    0    1]; %LITTLE_MANATEE
var_flag(7,:) = ...
          [    0     0    0    0     0    0    0    0    0    1]; %MANATEE
var_flag(8,:) = ...
          [    0     0    0    0     0    0    0    0    0    1]; %MYAKKA
var_flag(9,:) = ...
          [    0     0    0    0     0    0    0    0    0    1]; %PEACE
var_flag(10,:) = ...
          [    0     0    0    0     0    0    0    0    0    1]; %CALOOSAHATCHEE
var_flag(11,:) = ...
          [    0     0    0    0     0    0    0    0    0    1]; %LAKE

[calflow,caldat] = create_loadest_file(r, rname, varname(1:end-1), var_flag(:,1:end-1), n_river, out_date);

[calflow_sed,caldat_sed] = create_loadest_file_sed(r, rname, var_flag(:,end), n_river, out_date);

calflow(:,length(varname)) = calflow_sed;
caldat(:,length(varname)) = caldat_sed;

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
gen_bay_nutrients_ini('../Cruise_data',['./bay_initials_',num2str(year),'_09.xlsx']...
    ,grd,bay_mask,year,include_cruise,dep_range,time_range,time_range2,Ini_DON_method,bay_tss_flag);
end

%Step 7
%Read river sediment raw data
if(skip_step(7)==0)
    [tss_dat,tss_mean] = pre_preocess_river_tss('River_tss.xlsx',n_river);
    save(strcat('WA_',num2str(year),'.mat'),'tss_mean','tss_dat','-append');
end