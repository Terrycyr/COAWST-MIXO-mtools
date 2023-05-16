clear all;
grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');

r_i = [219 219 218 218 217 218 217 216 217 216 215 ...
    216 215 214 215 214 213 212 213 212 211 212 211 210 209 ...
    211 210 209 208 209 208 209 208];

r_j = [102 103 103 104 104 105 105 105 106 106 106 ...
    107 107 107 108 108 108 108 109 109 109 110 110 110 110 ...
    111 111 111 111 112 112 113 113];

for i=1:length(r_i)
    mask(r_i(i),r_j(i)) = 0;
end

figure;
contourf(mask');

save('removed_points.mat','r_i','r_j');

