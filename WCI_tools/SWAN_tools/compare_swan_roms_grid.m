clear all;

fn1 = 'SWAN_WFS_grid_OUTER.nc';
fn2 = 'C:\Users\cheny\Desktop\EcoHAB\Model_grid\ROMS_WFS_new.nc';

lon1 = ncread(fn1,'lon_rho');
lat1 = ncread(fn1,'lat_rho');
mask1 = ncread(fn1,'mask_rho');

lon2 = ncread(fn2,'lon_rho');
lat2 = ncread(fn2,'lat_rho');
mask2 = ncread(fn2,'mask_rho');

figure;
scatter(lon1(mask1>=0),lat1(mask1>=0),'r','filled');
hold on;
scatter(lon2(mask2>=0),lat2(mask2>=0),'g','filled');
