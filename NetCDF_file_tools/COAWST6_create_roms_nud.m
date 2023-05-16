clear all; close all;
grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
h = ncread(grd_name,'h');
N=gn.N;
fn = 'WFS_2021_nud.nc';

create_roms_netcdf_nudging_coef(fn,gn);

nud_2d = zeros(size(lon));
nud_3d = zeros(size(lon));

cff1_3d = 10;
cff1_2d = 10;
cff2 = 0;
n_gridpoints = 15;
decline_coef = 1.6;

for i=1:n_gridpoints
    nud0_2d(i) = (cff2+(n_gridpoints-i)^decline_coef*(cff1_2d-cff2)/(n_gridpoints-1)^decline_coef);
    nud0_3d(i) = (cff2+(n_gridpoints-i)^decline_coef*(cff1_3d-cff2)/(n_gridpoints-1)^decline_coef);
end

for k = 1:(n_gridpoints)

    nud_2d(1:end-k+1,k) = nud0_2d(k);
    nud_2d(end-k+1,k:end) = nud0_2d(k);  

    nud_3d(1:end-k+1,k) = nud0_3d(k);
    nud_3d(end-k+1,k:end) = nud0_3d(k);
end

for i=1:N
    nud3d(:,:,i) = nud_3d;
end

nud2d = nud_2d;

ncwrite(fn,'M3_NudgeCoef',nud3d);
ncwrite(fn,'M2_NudgeCoef',nud2d);
ncwrite(fn,'tracer_NudgeCoef',nud3d);

figure;
subplot(1,2,1);
nud_3d(mask==0) = NaN;
contourf(lon,lat,nud_3d,100,'linestyle','none');
subplot(1,2,2);
contourf(lon,lat,h,100,'linestyle','none');

figure;
plot(0:(length(nud0_2d)-1),nud0_2d,'k','LineWidth',2);refline(-cff1_2d/(n_gridpoints-1),cff1_2d)