clear all;close all;
fn =  '../Model_grid/ROMS_WFS_10river_grid.nc';
lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
lon_u = ncread(fn,'lon_u');
lat_u = ncread(fn,'lat_u');
lon_v = ncread(fn,'lon_v');
lat_v = ncread(fn,'lat_v');
angle = ncread(fn,'angle');
u = ncread('./WFS_2011_Uwind1.nc','Uwind');
v = ncread('./WFS_2011_Vwind1.nc','Vwind');
slp = ncread('./WFS_2011_Pair1.nc','Pair');
start_d = datenum(2011,1,1,0,0,0);

angle = -1*angle;
u2 = cos(angle).*u+sin(angle).*v;
v2 = -sin(angle).*u+cos(angle).*v;

for i=1:(62*4+1)
    figure(1)
    %contourf(lon,lat,slp(:,:,i));
    %hh = pcolor(lon,lat,sqrt(u(:,:,i).^2+v(:,:,i).^2));
    %set(hh,'linestyle','none');
    %contourf(lon,lat,v(:,:,i));
    quiver(lon(1:10:end,1:10:end),lat(1:10:end,1:10:end),u2(1:10:end,1:10:end,i),v2(1:10:end,1:10:end,i));
    title(datestr(start_d+i*(6/24)));
    colorbar;
    pause;
end



    
