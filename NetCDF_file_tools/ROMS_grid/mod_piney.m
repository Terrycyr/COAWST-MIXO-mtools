clear all; close all;
grd = 'ROMS_WFS_Piney.nc';
grd2 = 'ROMS_WFS_Piney_bay_mask.nc';

dep = ncread(grd,'h');
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
mask2 = ncread(grd2,'mask_rho');

dep_new = dep;
mask_new = mask;
mask_new2 = mask2;

mask_new(255,110) = 0;

mask_new(256,111) = 1;
dep_new(256,111) = 2;

mask_new2(255,110) = 0;
mask_new2(256,111) = 1;

rmask = mask_new;
[umask,vmask,pmask]=uvp_masks(rmask);

rmask2 = mask_new2;
[umask2,vmask2,pmask2]=uvp_masks(rmask2);

dep_new(mask==0) = 0.5;

ncwrite(grd,'h',dep_new);
ncwrite(grd,'mask_rho',rmask);
ncwrite(grd,'mask_u',umask);
ncwrite(grd,'mask_v',vmask);

ncwrite(grd2,'mask_rho',rmask2);
ncwrite(grd2,'mask_u',umask2);
ncwrite(grd2,'mask_v',vmask2);

figure;
mask_new(255,111) = 2;
pcolor(lon,lat,mask_new);
axis image;
