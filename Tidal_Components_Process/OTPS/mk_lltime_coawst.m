clear all;
grdname = '../../Model_grid/ROMS_WFS_10river_grid_v10.nc';
lon_whole=ncread(grdname,'lon_rho');
lat_whole=ncread(grdname,'lat_rho');
lon_whole_u=ncread(grdname,'lon_u');
lat_whole_u=ncread(grdname,'lat_u');
lon_whole_v=ncread(grdname,'lon_v');
lat_whole_v=ncread(grdname,'lat_v');

dt=1/24;
d1=datenum(2021,1,1,0,0,0);
d2=datenum(2022,1,1,0,0,0);

%South boundary
lon=lon_whole(:,1);
lat=lat_whole(:,1);
fid=fopen('lat_lon_time_s_el','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_u(:,1);
lat=lat_whole_u(:,1);
fid=fopen('lat_lon_time_s_u','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_v(:,1);
lat=lat_whole_v(:,1);
fid=fopen('lat_lon_time_s_v','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

%North boundary
lon=lon_whole(:,end);
lat=lat_whole(:,end);
fid=fopen('lat_lon_time_n_el','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_u(:,end);
lat=lat_whole_u(:,end);
fid=fopen('lat_lon_time_n_u','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_v(:,end);
lat=lat_whole_v(:,end);
fid=fopen('lat_lon_time_n_v','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);



%West boundary
lon=lon_whole(1,:);
lat=lat_whole(1,:);
fid=fopen('lat_lon_time_w_el','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_u(1,:);
lat=lat_whole_u(1,:);
fid=fopen('lat_lon_time_w_u','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_v(1,:);
lat=lat_whole_v(1,:);
fid=fopen('lat_lon_time_w_v','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

%East boundary
lon=lon_whole(end,:);
lat=lat_whole(end,:);
fid=fopen('lat_lon_time_e_el','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_u(end,:);
lat=lat_whole_u(end,:);
fid=fopen('lat_lon_time_e_u','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

lon=lon_whole_v(end,:);
lat=lat_whole_v(end,:);
fid=fopen('lat_lon_time_e_v','w');
for i=1:length(lat)
    for t=d1:dt:d2
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i),lon(i),datevec(t));
    end
end
fclose(fid);

%Initial field
lon=lon_whole;
lat=lat_whole;
fid=fopen('lat_lon_time_ini_el','w');
[r,c] = size(lon);
for i=1:r
    for j=1:c
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i,j),lon(i,j),datevec(d1));
    end
end
fclose(fid);

lon=lon_whole_u;%
lat=lat_whole_u;%
[r,c] = size(lon);
fid=fopen('lat_lon_time_ini_u','w');
for i=1:r
    for j=1:c
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i,j),lon(i,j),datevec(d1));
    end
end
fclose(fid);

lon=lon_whole_v;%
lat=lat_whole_v;%
[r,c] = size(lon);
fid=fopen('lat_lon_time_ini_v','w');
for i=1:r
    for j=1:c
        fprintf(fid,'%10.4f %10.4f %8d %4d %4d %4d %4d %4d\n',lat(i,j),lon(i,j),datevec(d1));
    end
end
fclose(fid);