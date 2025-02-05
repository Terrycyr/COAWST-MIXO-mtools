clear all;close all;
grd = './bio1dtest_grd.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = 21;
fn1 = 'biotest_Uwind1_bio.nc';
fn2 = 'biotest_Vwind1_bio.nc';
fn3 = 'biotest_swrad1_bio.nc';
fn4 = 'biotest_PAR_bio.nc';

[r,c] = size(lon);

time = 0:365;


%UWIND
Uwind = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn1,'Uwind');

%UWIND
Vwind = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn2,'Vwind');

%SSR
swrad = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn3,'swrad');

%PAR
PAR = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn4,'PAR');
