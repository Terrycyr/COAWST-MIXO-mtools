clear all;close all;

addpath(path,'C:\Users\cheny\Desktop\EcoHAB\NC_file_generation');

grd = './bio1dtest_grd.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = 21;
fn1 = 'biotest_Uwind1.nc';
fn2 = 'biotest_Vwind1.nc';
fn3 = 'biotest_Pair1.nc';
fn4 = 'biotest_Tair1.nc';
fn5 = 'biotest_Qair1.nc';
fn6 = 'biotest_swrad1.nc';
fn7 = 'biotest_lwrad1.nc';

[r,c] = size(lon);

time = 0:365;


%UWIND
Uwind = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn1,'Uwind');

%VWIND
Vwind = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn2,'Vwind');

%Pair
Pair = repmat(ones(size(lon)),1,1,length(time))*1025;
create_roms_forcings(lon,lat,time,fn3,'Pair');

%Tair
Tair = repmat(ones(size(lon)),1,1,length(time))*25;
create_roms_forcings(lon,lat,time,fn4,'Tair');

%Qair
Qair = repmat(ones(size(lon)),1,1,length(time))*1;
create_roms_forcings(lon,lat,time,fn5,'Qair');

%SSR
swrad = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn6,'swrad');

%SLR
lwrad = repmat(zeros(size(lon)),1,1,length(time));
create_roms_forcings(lon,lat,time,fn7,'lwrad');



