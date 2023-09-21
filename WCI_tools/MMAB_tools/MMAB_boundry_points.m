clear all
setup_nctoolbox
%lat = load('/home/terry/SWAN/grid/grid/PRE_lat.dat');
%lon = load('/home/terry/SWAN/grid/grid/PRE_lon.dat');
%dep = load('/home/terry/SWAN/grid/PRE_depth.dat');
lat = load('../COAWST/pre_coawstswan_2017/PRE_2017_lat.dat');
lon = load('../COAWST/pre_coawstswan_2017/PRE_2017_lon.dat');
dep = load('../COAWST/pre_coawstswan_2017/PRE_2017_depth.dat');
year = 2017;
month = 8;
R = 6371.004;
if(month<10)
    mydate = [num2str(year),'0',num2str(month)];
else
    mydate = [num2str(year),num2str(month)];
end
nc_dp = ncgeodataset(['./',num2str(year),'/multi_1.glo_30m.dp.',mydate,'.grb2']);
lat_dp = double(nc_dp{'lat'}(:));
lon_dp = double(nc_dp{'lon'}(:));
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
%    LatA = lat_w.*(pi/180);
%    LatB = lat(end,i).*(pi/180);
%    MLonA = lon_w.*(pi/180);
%    MLonB = lon(end,i).*(pi/180);
%    C = sin(LatA).*sin(LatB) + cos(LatA).*cos(LatB).*cos(MLonA-MLonB);
%    ds = R.*acos(C);
%    if(min(min(ds))<25)
%        if(m==0)
%            m=m+1;
%            [nb_y(m),nb_x(m)] = find(ds==min(min(ds)));
%        else
%            [rw,cw] = find(ds==min(min(ds)));
%            if(isempty(find(nb_y(:)==rw))||isempty(find(nb_x(:)==cw)))
%                m=m+1;
%                [nb_y(m),nb_x(m)] = find(ds==min(min(ds)));
%            end
%        end
%    end
end

delete('./MMAB_boudary_point.dat','wt+');
fid = fopen('./MMAB_boudary_point.dat','wt+');
l_w = length(wb_x);
l_e = length(eb_x);
l_s = length(sb_x);       
%l_n = length(nb_x);
for i=1:l_w
    fprintf(fid,'%5d%5d%6.1f%6.1f\n',wb_x(i),wb_y(i),w_lon(i),w_lat(i));
end

for i=1:l_e
    fprintf(fid,'%5d%5d%6.1f%6.1f\n',eb_x(i),eb_y(i),e_lon(i),e_lat(i));
end

for i=1:l_s
    fprintf(fid,'%5d%5d%6.1f%6.1f\n',sb_x(i),sb_y(i),s_lon(i),s_lat(i));
end

fclose(fid);

figure;
contourf(lon,lat,dep);
hold on;
scatter(w_lon,w_lat,'r');
hold on;
scatter(e_lon,e_lat,'g');
scatter(s_lon,s_lat,'b');