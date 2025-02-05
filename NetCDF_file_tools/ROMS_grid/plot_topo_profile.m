clear all;
grd = 'ROMS_WFS_10river_grid_v6.nc';
dep = ncread(grd,'h');
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');

[r,c] = size(mask);
for i=1:120
    x(i) = distance(lat(200,120),lon(200,120),lat(200,i),lon(200,i));
    y(i) = -1*mean(dep(200,i:i+1));
end

figure;
plot(x,y,'linewidth',2);
for i=1:21
    hold on;
    plot(x,y-y*(i/21),'color',[0.6 0.6 0.6],'linewidth',1);
end
set(gcf,'color','w');
axis([0 2.4 -3000 0]);
xlabel('Distance from coast (km)');


