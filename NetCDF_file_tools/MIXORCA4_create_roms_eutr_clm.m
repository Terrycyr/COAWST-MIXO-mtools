clear all; close all;
grd_name =  '../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
nfiles = 6;

year = 2005;
init_file = './Spinup/WFS_2005_ini_bio_mixo_sp.nc';

time_ref = datenum(year,6,1,0,0,0)-datenum(year,1,1);

date_clm = [0 732]-time_ref;
t_clim = length(date_clm);

%NO23
fn = 'WFS_2005_2006_clm_NO23_mixo.nc'; delete(fn);
create_roms_netcdf_clm_eutr(fn,gn,t_clim,'NO23');

no23_ini = ncread(init_file,'NO23');
no23_clm = repmat(no23_ini,1,1,1,2);

ncwrite(fn,'ocean_time',date_clm);
ncwrite(fn,'NO23_time',date_clm);
ncwrite(fn,'lon_rho',lon);
ncwrite(fn,'lat_rho',lat);
ncwrite(fn,'NO23',no23_clm);

%NH4T
fn = 'WFS_2005_2006_clm_NH4T_mixo.nc'; delete(fn);
create_roms_netcdf_clm_eutr(fn,gn,t_clim,'NH4T');

nh4t_ini = ncread(init_file,'NH4T');
nh4t_clm = repmat(nh4t_ini,1,1,1,2);

ncwrite(fn,'ocean_time',date_clm);
ncwrite(fn,'NH4T_time',date_clm);
ncwrite(fn,'lon_rho',lon);
ncwrite(fn,'lat_rho',lat);
ncwrite(fn,'NH4T',nh4t_clm);

%PO4T
fn = 'WFS_2005_2006_clm_PO4T_mixo.nc'; delete(fn);
create_roms_netcdf_clm_eutr(fn,gn,t_clim,'PO4T');

po4t_ini = ncread(init_file,'PO4T');
po4t_clm = repmat(po4t_ini,1,1,1,2);

ncwrite(fn,'ocean_time',date_clm);
ncwrite(fn,'PO4T_time',date_clm);
ncwrite(fn,'lon_rho',lon);
ncwrite(fn,'lat_rho',lat);
ncwrite(fn,'PO4T',po4t_clm);

%SIT
fn = 'WFS_2005_2006_clm_SIT_mixo.nc'; delete(fn);
create_roms_netcdf_clm_eutr(fn,gn,t_clim,'SIT');

sit_ini = ncread(init_file,'SIT');
sit_clm = repmat(sit_ini,1,1,1,2);

ncwrite(fn,'ocean_time',date_clm);
ncwrite(fn,'SIT_time',date_clm);
ncwrite(fn,'lon_rho',lon);
ncwrite(fn,'lat_rho',lat);
ncwrite(fn,'SIT',sit_clm);

%DO
fn = 'WFS_2005_2006_clm_DO_mixo.nc'; delete(fn);
create_roms_netcdf_clm_eutr(fn,gn,t_clim,'DO');

do_ini = ncread(init_file,'DO');
do_clm = repmat(do_ini,1,1,1,2);

ncwrite(fn,'ocean_time',date_clm);
ncwrite(fn,'DO_time',date_clm);
ncwrite(fn,'lon_rho',lon);
ncwrite(fn,'lat_rho',lat);
ncwrite(fn,'DO',do_clm);


