clear all; close all;
grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
h = ncread(grd_name,'h');
N=gn.N;
fn = 'WFS_2001_nud_bio.nc';

create_roms_netcdf_nudging_coef(fn,gn);

nud_2d = zeros(size(lon));
nud_3d = zeros(size(lon));

cff1_3d = 24;
cff1_2d = 24;
cff2 = 0;
n_gridpoints = 6;
s_gridpoints = 20;
e_gridpoints = 6;
w_gridpoints = 0;

for i=1:n_gridpoints
    nud0_n2d(i) = (cff2+((n_gridpoints-i)*(cff1_2d-cff2)/(n_gridpoints-1)));
    nud0_n3d(i) = (cff2+((n_gridpoints-i)*(cff1_3d-cff2)/(n_gridpoints-1)));
end


for k = 1:(n_gridpoints)
    nud_2d(k:end-k+1,end-k+1) = nud0_n2d(k);
    nud_3d(k:end-k+1,end-k+1) = nud0_n3d(k); 
end

% for k = 1:(n_gridpoints)
%     nud_2d(k,k:end-k+1) = nud0_2d(k);
%     nud_2d(k:end-k+1,k) = nud0_2d(k);
%     nud_2d(end-k+1,k:end-k+1) = nud0_2d(k);
%     nud_2d(k:end-k+1,end-k+1) = nud0_2d(k);  
% 
%     nud_3d(k,k:end-k+1) = nud0_3d(k);
%     nud_3d(k:end-k+1,k) = nud0_3d(k);
%     nud_3d(end-k+1,k:end-k+1) = nud0_3d(k);
%     nud_3d(k:end-k+1,end-k+1) = nud0_3d(k); 
% end

% for k = 1:(n_gridpoints)
%     nud_3d(:,k) = nud0_3d(k);
%     nud_2d(:,k) = nud0_2d(k);
% end

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