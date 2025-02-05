clear all; close all;

grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');

figure;
surf(lon,lat,-1*h,'linestyle','none');

ncwrite(grd_name,'Vstretching',1);
ncwrite(grd_name,'Vtransform',1);


