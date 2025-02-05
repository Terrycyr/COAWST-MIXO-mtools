clear all; close all;
grd = 'ROMS_WFS_10river_grid_v11.nc';
dep = ncread(grd,'h');
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');

dep_new = dep;
dep_new(dep_new<=3&mask==1) = 3;
dep_new(dep_new>500) = 500;

dep_new(mask==0) = 0.5;

mask_new = mask;

dep_new(158,159) = 4.5;
dep_new(158,160) = 3.5;
dep_new(159,159) = 3.5;

dep_new(164,148) = 3.5;
dep_new(163,148) = 3.5;
dep_new(181,141) = 3.5;
dep_new(185,137) = 3.5;
dep_new(182,139) = 3.5;
dep_new(187,131) = 3.7;
dep_new(187,132) = 3.7;
dep_new(186,133) = 3.7;
dep_new(186,135) = 3.7;
dep_new(186,134) = 3.7;
dep_new(190,127) = 3.7;

dep_new(133,160) = 3.5;
dep_new(133,161) = 4;
dep_new(133,162) = 3.5;

mask_new(241,105) = 0;
dep_new(241,105) = 0.5;

mask_new(234,99) = 1;
dep_new(234,99) = 3.2;

mask_new(284,97) = 1;
dep_new(284,97) = 4.5;

mask_new(285,97) = 1;
dep_new(285,97) = 4.5;


dep_new(283:285,98) = 3.8;

mask_new(350,103) = 0;
dep_new(350,103) = 0.5;

dep_new(351,103) = 10;

mask_new(352,103) = 1;
dep_new(352,103) = 8;

%dep_new(352,104:115) = dep_new(350,104:115);
%dep_new(353,104:115) = dep_new(349,104:115);
%dep_new(348,104:115) = dep_new(354,104:115);

mask_new(363,103) = 1;

dep_new(363,103) = 6;
dep_new(362,103) = 6;

dep_new(363,104) = 4.5;
dep_new(362,104) = 4.5;

dep_new(363,102) = 3.5;

mask_new(368,103) = 1;
dep_new(368,103) = 5.8;

mask_new(369,103) = 1;
dep_new(369,103) = 5.8;

dep_new(368,102) = 5;
dep_new(369,102) = 5;

dep_new(368,104) = 4.5;
dep_new(369,104) = 4.5;

dep_new(368,105) = 3.5;
dep_new(369,105) = 3.5;

dep_new(367,104) = 3.5;
dep_new(370,104) = 3.5;

dep_new(362,115) = 5.6;
dep_new(362,116) = 5.6773;
dep_new(361,115) = 5.6773;

mask_new(383,113:116) = 0;
dep_new(383,113:116) = 0.5;

dep_new(384,116) = 4;
dep_new(384:385,117) = 4;

dep_new(383,117:120) = 3.5;
dep_new(384,118) = 3.7;
dep_new(382,117:119) = 3.3;
dep_new(381,120) = 3.4;
dep_new(382,120) = 4;
dep_new(382,121) = 5;
dep_new(382,122) = 4;
dep_new(383,122) = 4;

dep_new(384,120) = 3.53;
dep_new(385,118) = 3.9285;

tmp = dep_new(396:end,:);
tmp(tmp<4) = 4;
dep_new(396:end,:)=tmp;

dep_new(mask==0) = 0.5;

rmask = mask_new;
[umask,vmask,pmask]=uvp_masks(rmask);

ncwrite(grd,'h',dep_new);
ncwrite(grd,'mask_rho',rmask);
ncwrite(grd,'mask_u',umask);
ncwrite(grd,'mask_v',vmask);

figure;
