clear all;close all;

addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

grd_name =  'ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
bay_mask = ncread('ROMS_WFS_new_bay_mask.nc','mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
N = length(ncread(grd_name,'Cs_r'));

dep(bay_mask==1) = NaN;

fid = fopen('floats_upwelling.in','wt+');
fid0 = fopen('mod_floats_upwelling.in','rt+');

year = 2022;

ini_date = datenum(year,1,1);

release_time = [datenum(2022,9,27),datenum(2022,9,30),datenum(2022,10,7)];
release_rounds = 3;

sparse_coef = 4;
sparse_flag = ones(size(dep));
sparse_flag(1:sparse_coef:end,1:sparse_coef:end) = 0;

pos_flag = dep>10&dep<45;
pos_flag(250,92,:) = 0;
pos_flag(251,93,:) = 0;
pos_flag(1:190,:) = 0.;
pos_flag(350:end,:) = 0.;
pos_flag(sparse_flag==1) = 0;
tmp1 = lon(pos_flag);
tmp2 = lat(pos_flag);

release_pos{1} = [tmp1';tmp2'];
release_pos{2} = release_pos{1};
release_pos{3} = release_pos{1};
release_pos{4} = release_pos{1};

obs_pos = [-83.2979  -83.4104  -83.2474; 27.7014   27.6497   27.6497];

for rd = 1:release_rounds
    for i = 1:size(release_pos{rd},2)
        [release_ij{rd}(1,i),release_ij{rd}(2,i)] = find_ij(release_pos{rd}(1,i),release_pos{rd}(2,i),lon,lat,mask);

    end
end

release_layer = [1 11 21];

figure;
contourf(lon,lat,dep,'LineStyle','none');
hold on;
scatter(release_pos{1}(1,:),release_pos{1}(2,:),'black','filled');
hold on;
scatter(obs_pos(1,:),obs_pos(2,:),'r','filled');


%Nested grid number [ng]. This must be 1 in non-nested applications.
G=1;

% Initial horizontal location (Fx0 and Fy0) coordinate type.
% If C = 1, then Fx0 is longitude (west values are negative) 
%and Fx0 is latitude (south values are negative).
C=0;

% If T = 1, float(s) will be 3D Lagrangrian particles. 
% If T = 2, float(s) will be isobaric particles . 
% If T = 3, float(s) will be geopotential (constant depth) particles.
T=1;

% Number of floats to be released at the (Fx0,Fy0,Fz0) location
N=1;

N_total = 0;
for rd = 1:release_rounds
    N_total = N_total+size(release_pos{rd},2).*length(release_layer);
end


while(~feof(fid0))
    tmp = fgetl(fid0);
    tmp2 = strtrim(tmp);
    if(length(tmp2)>7)
        if(strcmp(tmp2(1:7),'NFLOATS'))
            fprintf(fid,'     NFLOATS == %d',N_total);
        elseif(strcmp(tmp2(1:3),'POS'))
            fprintf(fid,'%s\n',tmp);
            fprintf(fid,'\n');
            for rd = 1:release_rounds
                Ft0 = release_time(rd)-ini_date;
                Fx0 = release_ij{rd}(1,:);
                Fy0 = release_ij{rd}(2,:);
                Fdt = 0;
                Fdx = 0;
                Fdy = 0;
                Fdz = 0;
                for flt = 1:length(Fx0)
                    for ly = 1:length(release_layer)
                        Fz0 = release_layer(ly);
                        fprintf(fid,'%7d%3d%3d%3d%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0%6.1fd0\n'...
                            ,G,C,T,N,Ft0,Fx0(flt),Fy0(flt),Fz0,Fdt,Fdx,Fdy,Fdz);
                    end
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






