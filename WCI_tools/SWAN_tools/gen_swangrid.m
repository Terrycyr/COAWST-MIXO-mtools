clear all;
grd = 'ROMS_WFS_new.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
dep = ncread(grd,'h');
mask = ncread(grd,'mask_rho');
lon = lon';
lat = lat';
dep = dep';
mask = mask';
dep(mask==0)=-999.999;

fid = fopen('./WFS_swan_grd.dat','wt+');
fid2 = fopen('./WFS_swan_depth.dat','wt+');
fid3 = fopen('./WFS_swan_lon.dat','wt+');
fid4 = fopen('./WFS_swan_lat.dat','wt+');
fid5 = fopen('./WFS_swan_coord.dat','wt+');
[r,c] = size(lat);

for i=1:r
    fprintf(fid,'%13.6f',lon(i,:));
    fprintf(fid3,'%13.6f',lon(i,:));
    fprintf(fid2,'%10.3f',dep(i,:));
    fprintf(fid,'\n');
    fprintf(fid2,'\n');
    fprintf(fid3,'\n');
    for j=1:length(lon(i,:))
        fprintf(fid5,'%13.6f',lon(i,j));
        fprintf(fid5,'\n');
    end
end

for i=1:r
    fprintf(fid,'%13.6f',lat(i,:));
    fprintf(fid4,'%13.6f',lat(i,:));
    fprintf(fid,'\n');
    fprintf(fid4,'\n');
    for j=1:length(lat(i,:))
        fprintf(fid5,'%13.6f',lat(i,j));
        fprintf(fid5,'\n');
    end
end
fclose(fid);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
