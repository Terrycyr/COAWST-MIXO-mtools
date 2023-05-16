clear all;
addpath C:\Users\cheny\Desktop\EcoHAB\self_functions
grd = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd,'Cs_r'));

fn1 = 'WFS_2010_swrad1.nc';

load('C:\Users\cheny\Desktop\EcoHAB\PAR\erdMGpar01day_c46c_e27d_dbdf.mat')

time = erdMGpar01day.time;
tvec = datevec(time/3600/24+datenum(1970,1,1));
par = mean(mean(double(erdMGpar01day.par),4),3);

pos = find(tvec(:,1)==2010);
par_lon = mean(erdMGpar01day.longitude);
par_lat = mean(erdMGpar01day.latitude);

[i,j] = find_ij(par_lon,par_lat,lon,lat,mask);

swrad = squeeze(double(ncread(fn1,'swrad',[i,j,1],[1,1,inf])));
for i=1:365
    swrad_m(i) = mean(swrad([1:4]+4*(i-1)));
end

par2 = par(pos)*1e6/24/3600; %Einsteins m-2 d-1 -> µmoles m-2 s-1
figure;
plot(1:365,par2,linspace(1,365,length(swrad)),swrad*1.903+89.12)


