%Modified by terry, 2021-11-07
clear all;

grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));

out_path = strcat('./RCA/');
if(~exist(out_path))
    mkdir(out_path);
end

fn = [out_path,'WFS_eutr_2001_bry.nc'];

delete(fn);

n_river = 11;
dir_flag = [1 1 0 1]; %S N W E

date_out = datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,24,0,0);
year = datevec(date_out(1));
year = year(1);
t = length(date_out);
create_rca_netcdf_bndry_eutr(fn,gn,t,n_river,dir_flag)
const_initialize(fn,0.);

load(strcat('../Water_Atlas/','WA_river_bnd_',num2str(year),'.mat'));
load(strcat('../GOM_preprocessing/','ocean_nutrients_bnd_',num2str(year),'.mat'));

%
SAL_south = s_sal;
DO_south = s_do;
NO23_south = s_no3/10;
NH4T_south = s_nh4;
PO4T_south = s_po4;
SIT_south = s_si;

SAL_north = n_sal;
DO_north = n_do;
NO23_north = n_no3/10;
NH4T_north = n_nh4;
PO4T_north = n_po4;
SIT_north = n_si;

SAL_east = e_sal;
DO_east = e_do;
NO23_east = e_no3*0.;
NH4T_east = e_nh4*0.;
PO4T_east = e_po4*0.;
SIT_east = e_si*0.;

river_NO23(1:3,:,:) = 0.1*river_NO23(1:3,:,:);
river_NH4T(1:3,:,:) = 0.1*river_NH4T(1:3,:,:);
river_PO4T(1:3,:,:) = 0.1*river_PO4T(1:3,:,:);

ncwrite(fn,'bry_time',date_out-date_out(1));
ncwrite(fn,'IBCOPT',2);
ncwrite(fn,'IBCPWLOPT',1);
ncwrite(fn,'bry_time',date_out-date_out(1));
%
ncwrite(fn,'river_LDOC',river_LDOC);
ncwrite(fn,'river_LDOC',river_LDOC);
ncwrite(fn,'river_LDON',river_LDON);
ncwrite(fn,'river_LDOP',river_LDOP);
ncwrite(fn,'river_RDOC',river_RDOC);
ncwrite(fn,'river_RDON',river_RDON);
ncwrite(fn,'river_RDOP',river_RDOP);
ncwrite(fn,'river_NH4T',river_NH4T);
ncwrite(fn,'river_NO23',river_NO23);
ncwrite(fn,'river_PO4T',river_PO4T);
ncwrite(fn,'river_LPOC',river_LPOC);
ncwrite(fn,'river_LPON',river_LPON);
ncwrite(fn,'river_LPOP',river_LPOP);
ncwrite(fn,'river_RPOC',river_RPOC);
ncwrite(fn,'river_RPON',river_RPON);
ncwrite(fn,'river_RPOP',river_RPOP);
ncwrite(fn,'river_SAL',river_SAL);
ncwrite(fn,'river_SIT',river_SIT);
ncwrite(fn,'river_DO',river_DO);

%south
ncwrite(fn,'SAL_south',SAL_south,[2 1 1]);
ncwrite(fn,'SIT_south',SIT_south,[2 1 1]);
ncwrite(fn,'DO_south',DO_south,[2 1 1]);
ncwrite(fn,'NO23_south',NO23_south,[2 1 1]);
ncwrite(fn,'NH4T_south',NH4T_south,[2 1 1]);
ncwrite(fn,'PO4T_south',PO4T_south,[2 1 1]);

% %north
% ncwrite(fn,'SAL_north',SAL_north,[2 1 1]);
% ncwrite(fn,'SIT_north',SIT_north,[2 1 1]);
% ncwrite(fn,'DO_north',DO_north,[2 1 1]);
% ncwrite(fn,'NO23_north',NO23_north,[2 1 1]);
% ncwrite(fn,'NH4T_north',NH4T_north,[2 1 1]);
% ncwrite(fn,'PO4T_north',PO4T_north,[2 1 1]);

%east
ncwrite(fn,'SAL_east',SAL_east,[2 1 1]);
ncwrite(fn,'SIT_east',SIT_east,[2 1 1]);
ncwrite(fn,'DO_east',DO_east,[2 1 1]);
ncwrite(fn,'NO23_east',NO23_east,[2 1 1]);
ncwrite(fn,'NH4T_east',NH4T_east,[2 1 1]);
ncwrite(fn,'PO4T_east',PO4T_east,[2 1 1]);