clear all; close all;
grd = 'ROMS_WFS_10river_grid_v3.nc';
dep = ncread(grd,'h');
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');

dep(mask==0) = NaN;

fn = '../HYCOM/2010/HYCOM_20100101.nc';
hycom_lon = double(ncread(fn,'lon'));
hycom_lat = double(ncread(fn,'lat'));
hycom_dep = double(ncread(fn,'depth'));
u  =double(ncread(fn, 'u'));

[a,b,c,d] = size(u);
dv = ~isnan(u);
dv2 = squeeze(sum(dv,3));

dv2(dv2==0) = 1;

for i = 1:a
    for j =1:b
        topo(i,j) = hycom_dep(dv2(i,j,1));
    end
end

topo2 = griddata(hycom_lon,hycom_lat,topo,lon,lat);

nud = ones(size(topo2));
nud0 = linspace(0,1,6);

[a,b,c,d] = size(topo2);



for k = 1:5
    nud(k,k:end-k+1) = nud0(k);
    nud(k:end-k+1,k) = nud0(k);
    nud(end-k+1,k:end-k+1) = nud0(k);
    nud(k:end-k+1,end-k+1) = nud0(k);  
end

dep_new = dep.*nud+topo2.*(1-nud);
figure;
contourf(lon,lat,dep_new);

%ncwrite(grd,'h',dep_new);
