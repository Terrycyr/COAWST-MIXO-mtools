clear all;
lon = ncread('./pre_2017_grd.nc','lon_rho');
lat = ncread('./pre_2017_grd.nc','lat_rho');
x = ncread('./pre_2017_grd.nc','x_rho');
y = ncread('./pre_2017_grd.nc','y_rho');
dep = ncread('./pre_2017_grd.nc','h');
mask = ncread('./pre_2017_grd.nc','mask_rho');
lon = lon';
lat = lat';
x = x';
y = y';
dep = dep';
mask = mask';
dep(mask==0)=-999.999;

fid = fopen('./PRE_2017_grd.dat','wt+');
fid2 = fopen('./PRE_2017_depth.dat','wt+');
fid3 = fopen('./PRE_2017_lon.dat','wt+');
fid4 = fopen('./PRE_2017_lat.dat','wt+');
fid5 = fopen('./PRE_2017_xy.dat','wt+');
[r,c] = size(lat);

for i=1:r
    fprintf(fid,'%13.6f',lon(i,:));
    fprintf(fid3,'%13.6f',lon(i,:));
    fprintf(fid5,'%12.1f',x(i,:));
    fprintf(fid2,'%10.3f',dep(i,:));
    fprintf(fid,'\n');
    fprintf(fid2,'\n');
    fprintf(fid3,'\n');
    fprintf(fid5,'\n');
end

for i=1:r
    fprintf(fid,'%13.6f',lat(i,:));
    fprintf(fid4,'%13.6f',lat(i,:));
    fprintf(fid5,'%12.1f',y(i,:));
    fprintf(fid,'\n');
    fprintf(fid4,'\n');
    fprintf(fid5,'\n');
end
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
