clear all; close all;
grd_name =  '../../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
h = ncread(grd_name,'h');
N=gn.N;
fn = 'WFS_2005_nud_bio_mixo_sp.nc';

create_roms_netcdf_nudging_coef(fn,gn);

nud_2d = zeros(size(lon));
nud_3d = zeros(size(lon));

cff1_3d = 4;
cff1_2d = 4;
cff2 = 1/30;
n_gridpoints = 20;

for i=1:n_gridpoints
    nud0_2d(i) = (cff2+(n_gridpoints-i)^2*(cff1_2d-cff2)/(n_gridpoints-1)^2);
    nud0_3d(i) = (cff2+(n_gridpoints-i)^2*(cff1_3d-cff2)/(n_gridpoints-1)^2);
end

for k = 1:(n_gridpoints)

    nud_2d(1:end-k+1,k) = nud0_2d(k);
    nud_3d(1:end-k+1,k) = nud0_3d(k);    
    
    if(k<10)
    nud_2d(end-k+1,k:end) = nud0_2d(k);  
    nud_3d(end-k+1,k:end) = nud0_3d(k);
    end
end

for i=1:N
    nud3d(:,:,i) = nud_3d;
end

nud2d = nud_2d;

ncwrite(fn,'M3_NudgeCoef',nud3d);
ncwrite(fn,'M2_NudgeCoef',nud2d);
ncwrite(fn,'tracer_NudgeCoef',nud3d);

figure;
nud_3d(mask==0) = NaN;
hp = pcolor(lon,lat,nud_3d);
set(hp,'linestyle','none');