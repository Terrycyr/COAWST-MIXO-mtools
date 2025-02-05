clear all

lat = load('WFS_swan_outer_lat.dat');
lon = load('WFS_swan_outer_lon.dat');
dep = load('WFS_swan_outer_depth.dat');
year = 2022;
month = 1;
R = 6371.004;

%list for included WW3 files
hycom_dir = ['./',num2str(year),'/'];
fname = dir([hycom_dir 'WW3*.nc']);
for i = 1:size(fname,1)
    filename{i,1} = strcat(hycom_dir,fname(i).name);
end

nc_dp = filename{1};
lat_dp = double(ncread(nc_dp,'lat'));
lon_dp = double(ncread(nc_dp,'lon'));
[r,c] = size(lat);
[lon_w,lat_w] = meshgrid(lon_dp,lat_dp);

% W,E
n = 0;
m = 0;
for i=1:r
    if(dep(i,1)~=-999.999)
        LatA = lat_w.*(pi/180);
        LatB = lat(i,1).*(pi/180);
        MLonA = lon_w.*(pi/180);
        MLonB = lon(i,1).*(pi/180);
        C = sin(LatA).*sin(LatB) + cos(LatA).*cos(LatB).*cos(MLonA-MLonB);
        ds = R.*acos(C);
        if(min(min(ds))<25)
            if(n==0)
                n=n+1;
                [wb_y(n),wb_x(n)] = find(ds==min(min(ds)));
                w_lat(n) = lat_w(wb_y(n),wb_x(n));
                w_lon(n) = lon_w(wb_y(n),wb_x(n));
            else
                [rw,cw] = find(ds==min(min(ds)));
                if(isempty(find(wb_y(:)==rw))||isempty(find(wb_x(:)==cw)))
                    n=n+1;
                    [wb_y(n),wb_x(n)] = find(ds==min(min(ds)));
                    w_lat(n) = lat_w(wb_y(n),wb_x(n));
                    w_lon(n) = lon_w(wb_y(n),wb_x(n));
                end
            end
        end
    end
    if(dep(i,end)~=-999.999)
        LatA = lat_w.*(pi/180);
        LatB = lat(i,end).*(pi/180);
        MLonA = lon_w.*(pi/180);
        MLonB = lon(i,end).*(pi/180);
        C = sin(LatA).*sin(LatB) + cos(LatA).*cos(LatB).*cos(MLonA-MLonB);
        ds = R.*acos(C);
        if(min(min(ds))<25)
            if(m==0)
                m=m+1;
                [eb_y(m),eb_x(m)] = find(ds==min(min(ds)));
                e_lat(m) = lat_w(eb_y(m),eb_x(m));
                e_lon(m) = lon_w(eb_y(m),eb_x(m));
            else
                [rw,cw] = find(ds==min(min(ds)));
                if(isempty(find(eb_y(:)==rw))||isempty(find(eb_x(:)==cw)))
                    m=m+1;
                    [eb_y(m),eb_x(m)] = find(ds==min(min(ds)));
                    e_lat(m) = lat_w(eb_y(m),eb_x(m));
                    e_lon(m) = lon_w(eb_y(m),eb_x(m));
                end
            end
        end
    end    
end

% N,S
n = 0;
m = 0;
for i=1:c
    if(dep(1,i)~=-999.999)
        LatA = lat_w.*(pi/180);
        LatB = lat(1,i).*(pi/180);
        MLonA = lon_w.*(pi/180);
        MLonB = lon(1,i).*(pi/180);
        C = sin(LatA).*sin(LatB) + cos(LatA).*cos(LatB).*cos(MLonA-MLonB);
        ds = R.*acos(C);
        if(min(min(ds))<25)
            if(n==0)
                n=n+1;
                [sb_y(n),sb_x(n)] = find(ds==min(min(ds)));
                s_lat(n) = lat_w(sb_y(n),sb_x(n));
                s_lon(n) = lon_w(sb_y(n),sb_x(n));
            else
                [rw,cw] = find(ds==min(min(ds)));
                if(isempty(find(sb_y(:)==rw))||isempty(find(sb_x(:)==cw)))
                    n=n+1;
                    [sb_y(n),sb_x(n)] = find(ds==min(min(ds)));
                    s_lat(n) = lat_w(sb_y(n),sb_x(n));
                    s_lon(n) = lon_w(sb_y(n),sb_x(n));
                end
            end
        end
    end
    if(dep(end,i)~=-999.999)
        LatA = lat_w.*(pi/180);
        LatB = lat(end,i).*(pi/180);
        MLonA = lon_w.*(pi/180);
        MLonB = lon(end,i).*(pi/180);
        C = sin(LatA).*sin(LatB) + cos(LatA).*cos(LatB).*cos(MLonA-MLonB);
        ds = R.*acos(C);
        if(min(min(ds))<25)
            if(m==0)
                m=m+1;
                [nb_y(m),nb_x(m)] = find(ds==min(min(ds)));
                n_lat(m) = lat_w(nb_y(m),nb_x(m));
                n_lon(m) = lon_w(nb_y(m),nb_x(m));
            else
                [rw,cw] = find(ds==min(min(ds)));
                if(isempty(find(nb_y(:)==rw))||isempty(find(nb_x(:)==cw)))
                    m=m+1;
                    [nb_y(m),nb_x(m)] = find(ds==min(min(ds)));
                    n_lat(m) = lat_w(nb_y(m),nb_x(m));
                    n_lon(m) = lon_w(nb_y(m),nb_x(m));
                end
            end
        end
    end
end

delete('./WW3_boudary_point.dat','wt+');
fid = fopen('./WW3_boudary_point.dat','wt+');
if(exist("wb_x","var"))
l_w = length(wb_x);
else
    l_w = [];
end

if(exist("eb_x","var"))
l_e = length(eb_x);
else
    l_e = [];
end

if(exist("sb_x","var"))
l_s = length(sb_x);
else
    l_s = [];
end

if(exist("nb_x","var"))
l_n = length(nb_x);
else
    l_n = [];
end

for i=1:l_w
    fprintf(fid,'%5d%5d%6.1f%6.1f\n',wb_x(i),wb_y(i),w_lon(i),w_lat(i));
end

for i=1:l_e
    fprintf(fid,'%5d%5d%6.1f%6.1f\n',eb_x(i),eb_y(i),e_lon(i),e_lat(i));
end

for i=1:l_s
    fprintf(fid,'%5d%5d%6.1f%6.1f\n',sb_x(i),sb_y(i),s_lon(i),s_lat(i));
end

for i=1:l_n
    fprintf(fid,'%5d%5d%6.1f%6.1f\n',nb_x(i),nb_y(i),n_lon(i),n_lat(i));
end

fclose(fid);

figure;
contourf(lon,lat,dep);
hold on;
if(exist("wb_x","var"))
    scatter(w_lon,w_lat,'r');
end
hold on;
if(exist("eb_x","var"))
    scatter(e_lon,e_lat,'g');
end
hold on;
if(exist("sb_x","var"))
    scatter(s_lon,s_lat,'b');
end
hold on;
if(exist("nb_x","var"))
    scatter(n_lon,n_lat,'y');
end