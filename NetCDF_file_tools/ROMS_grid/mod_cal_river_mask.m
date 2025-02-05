clear all; close all;
grd = 'ROMS_WFS_new.nc';
grd2 = 'ROMS_WFS_new_mask.nc';

dep = ncread(grd,'h');
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
mask2 = ncread(grd2,'mask_rho');

dep_new = dep;
mask_new = mask;
mask_new2 = mask2;

mask_new(376,119) = 1;
dep_new(376,119) = 3;

mask_new(377,116) = 1;
dep_new(377,116) = 3;

mask_new(378,116) = 1;
dep_new(378,116) = 3;

mask_new2(376,119) = 1;
mask_new2(377,116) = 1;
mask_new2(378,116) = 1;

rmask = mask_new;
[umask,vmask,pmask]=uvp_masks(rmask);

rmask2 = mask_new2;
[umask2,vmask2,pmask2]=uvp_masks(rmask2);

ncwrite(grd,'h',dep_new);
ncwrite(grd,'mask_rho',rmask);
ncwrite(grd,'mask_u',umask);
ncwrite(grd,'mask_v',vmask);

ncwrite(grd2,'mask_rho',rmask2);
ncwrite(grd2,'mask_u',umask2);
ncwrite(grd2,'mask_v',vmask2);

figure;
pcolor(lon,lat,mask_new);
axis image;
