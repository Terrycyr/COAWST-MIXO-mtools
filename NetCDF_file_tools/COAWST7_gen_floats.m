clear all;

grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v10.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));

fid = fopen('floats_upwelling.in','wt+');
fid0 = fopen('mod_floats_upwelling.in','rt+');

year = 2002;

release_isobath = 75; %meter
release_depth = 1.5; %meter from bottom
release_month = 1:12;

search_range = 5; %meter
select_range = 20:360; % grid - i index

[r,c] = size(lon);
k=0;
for i=1:r
    tmp = dep(i,:);
    dis = abs(tmp-release_isobath);
    dis(isnan(dis)) = 10000;
    dis(dis>search_range) = 10000;
    if(length(find(dis==min(dis)))==1)
        k=k+1;
        s_j(k) = find(dis==min(dis));
        s_i(k) = i;
    end
end

figure;
contourf(lon,lat,dep,'LineStyle','none');
hold on;
for i=select_range
    scatter(lon(s_i(i),s_j(i)),lat(s_i(i),s_j(i)),'black','filled');
    hold on;
end

%Nested grid number [ng]. This must be 1 in non-nested applications.
G=1;

% Initial horizontal location (Fx0 and Fy0) coordinate type.
% If C = 0, then Fx0 is longitude (west values are negative) 
%and Fx0 is latitude (south values are negative).
C=0;

% If T = 1, float(s) will be 3D Lagrangrian particles. 
% If T = 2, float(s) will be isobaric particles . 
% If T = 3, float(s) will be geopotential (constant depth) particles.
T=1;

% Number of floats to be released at the (Fx0,Fy0,Fz0) location
N=1;

Fx0 = s_i(select_range);
Fy0 =  s_j(select_range);
Fz0 = release_depth;
Fdt = 0;
Fdx = 0;
Fdy = 0;
Fdz = 0;

while(~feof(fid0))
    tmp = fgetl(fid0);
    tmp2 = strtrim(tmp);
    if(length(tmp2)>7)
        if(strcmp(tmp2(1:7),'NFLOATS'))
            fprintf(fid,'     NFLOATS == %d',length(release_month)*length(Fx0));
        elseif(strcmp(tmp2(1:3),'POS'))
            fprintf(fid,'%s\n',tmp);
            fprintf(fid,'\n');
            for month = release_month
                Ft0 = datenum(year,month,1)-datenum(year,1,1);
                for flt = 1:length(Fx0)
                    fprintf(fid,'%7d%3d%3d%3d%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0\n'...
                        ,G,C,T,N,Ft0,Fx0(flt),Fy0(flt),Fz0,Fdt,Fdx,Fdy,Fdz);
                end
            end
        else
            fprintf(fid,'%s\n',tmp);
        end
    else
        fprintf(fid,'%s\n',tmp);
    end  
end

fclose(fid);
fclose(fid0);






