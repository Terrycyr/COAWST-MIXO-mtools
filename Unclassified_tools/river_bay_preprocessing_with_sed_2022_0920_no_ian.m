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
skip_step = [0 0 0 0 0];

hurricane_lanfall = datenum(2022,9,28,20,35,00);
replace_len_before = 2;
replace_len_after = 15;
%Make sure selected right file for initials, might be different for ROMS and RCA.

fn_WA = strcat('WA_',num2str(year),'_no_IAN.mat');
fn_calib = strcat('calib_',num2str(year),'_no_IAN.mat');
fn_loadest = strcat('loadest_',num2str(year),'_no_IAN.mat');
fn_rout = strcat('WA_river_bnd_',num2str(year),'_no_IAN.mat');

%Step 1
%Read river raw data
if(skip_step(1)==0)
    [nutri_dat,nutri_mean] = pre_process_river_nutrients('2022-no-ian.xlsx',varname(1:end-1),n_river,on_flag);  %!!!!
    save(fn_WA,'nutri_mean','nutri_dat');
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

%Remove Ian signial
for i=1:size(r,2)
    tmp = r(i).data;
    rep_pos = find((tmp(:,4)>=hurricane_lanfall-replace_len_before)...
    &(tmp(:,4)<=hurricane_lanfall+replace_len_after));
    tmp(rep_pos,:) = [];
    r(i).data = tmp;
end

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
        read_plot_loadest_result(rname, varname, var_flag, n_river...
        , out_date, './fig/',fn_WA,fn_calib);

    save(fn_loadest,'result_conc','result_flow');
end

%Step 5 river boundaries
if(skip_step(5)==0)
    gen_river_nutrients(grd, varname, var_flag, n_river, out_date,fn_WA,fn_loadest,fn_rout);
end