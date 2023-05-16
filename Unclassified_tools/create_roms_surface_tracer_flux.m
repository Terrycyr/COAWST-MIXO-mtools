clear all;close all;
grd = '../Model_grid/ROMS_WFS_10river_grid_v10.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd,'Cs_r'));
fn1 = 'WFS_2001_Dye1.nc';

origin_date = datenum(2001,1,1,0,0,0):6/24:datenum(2001,12,31,24,0,0);

out_date1 = datenum(2001,1,1,0,0,0):6/24:datenum(2001,12,31,24,0,0);
out_date2 = [];

time1 = out_date1-out_date1(1);
time2 = out_date2-out_date1(1);

[r,c]  = size(lon);

dye01 = zeros(r,c,length(time1));
dye_01_sflux = dye01;
if(sum(sum(sum(isnan(dye_01_sflux))))>0)
    pause;
end

create_roms_forcings(lon,lat,time1,fn1,'dye_01_sflux');
