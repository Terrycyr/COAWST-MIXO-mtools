clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check following variables when setting up the script
% 1. Grid name.
% 2. Year
% 3. nontidal_dataset_dir and nontidal_dataset_source
% 4. Output file name
% 5. sed_flag, dye_flag and river_bnd
% 6. adjust_tide and adjust_el
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\NC_file_generation');

grd_name =  '../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
NCS = 0;
NNS = 0;
dye_flag=0;
river_bnd = 0;
adjust_tide = 1.;
%adjust_el = mean(mean(s_el_final(:,1:l_period)));
adjust_el = 0.35;

year = 2003;
nontidal_dataset_dir = '../Non_tidal_component_preprocessing/HYCOM/';
%nontidal_dataset_source = 'GLBy0.08/expt_56.3_57.2';
nontidal_dataset_source = 'GLBv0.08/expt_53.X';

fn = ['WFS_',num2str(year),'_bry.nc'];
date_out3d = datenum(year,1,1,0,0,0):6/24:datenum(year,12,31,24,0,0);
date_out2d = datenum(year,1,1,0,0,0):1/24:datenum(year,12,31,24,0,0);


if(dye_flag==1)
    date_out_dye = datenum(year,1,1,0,0,0):1:datenum(year,12,31,24,0,0);
end
t_clim = length(date_out3d);
t_clim2 = length(date_out2d);
if(dye_flag==1)
    t_clim3 = length(date_out_dye);
else
    t_clim3 = [];
end

create_roms_netcdf_bndry_mwUL(fn,gn,t_clim,t_clim2,t_clim3,NCS,NNS,nontidal_dataset_source);

[r,c] = size(lon);
zeta_time = date_out2d - date_out2d(1);
v2d_time = date_out2d - date_out2d(1);
v3d_time = date_out3d - date_out3d(1);
salt_time = date_out3d - date_out3d(1);
temp_time = date_out3d - date_out3d(1);
if(dye_flag==1)
    dye_time = date_out_dye - date_out_dye(1);
end

if(NCS>0)
    mud_time = date_out3d - date_out3d(1);
end

if(NNS>0)
    sand_time = date_out3d - date_out3d(1);
end

load(strcat('../Tidal_component_preprocessing/tide_bnd_',num2str(year),'.mat'));
load(strcat(nontidal_dataset_dir,'non_tidal_bnd_',num2str(year),'.mat'));

if(dye_flag==1&&year~=2021)
    load(strcat('../GOM_preprocessing/ocean_nutrients_bnd_',num2str(year),'.mat'));
end

if(river_bnd==1)
    load(strcat('./river_raw/river_bnd_',num2str(year),'.mat'));
    %load('../Water_Atlas/WAtlas_river_bnd.mat');
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
for i=1:gn.N
    if(dye_flag==1&&year~=2021)
        s_tracer01(:,gn.N-i+1,:) = s_no3(:,i,:);
        n_tracer01(:,gn.N-i+1,:) = n_no3(:,i,:);
        e_tracer01(:,gn.N-i+1,:) = e_no3(:,i,:);
    end
% s_temp2(:,gn.N-i+1,:) = s_temp(:,i,:);
% n_temp2(:,gn.N-i+1,:) = n_temp(:,i,:);
% e_temp2(:,gn.N-i+1,:) = e_temp(:,i,:);
% s_s2(:,gn.N-i+1,:) = s_sal(:,i,:);
% n_s2(:,gn.N-i+1,:) = n_sal(:,i,:);
% e_s2(:,gn.N-i+1,:) = e_sal(:,i,:);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tide_n_el = tide_n_el*adjust_tide;
tide_s_el = tide_s_el*adjust_tide;
%tide_w_el = tide_w_el*adjust_tide;
tide_e_el = tide_e_el*adjust_tide;

tide_n_v = tide_n_v*adjust_tide;
tide_s_v = tide_s_v*adjust_tide;
%tide_w_v = tide_w_v*adjust_tide;
tide_e_v = tide_e_v*adjust_tide;

tide_n_u = tide_n_u*adjust_tide;
tide_s_u = tide_s_u*adjust_tide;
%tide_w_u = tide_w_u*adjust_tide;
tide_e_u = tide_e_u*adjust_tide;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n_el_final = tide_n_el+n_el_out;
s_el_final =  tide_s_el+s_el_out;
% w_el_final = tide_w_el+w_el_out;
 e_el_final =  tide_e_el+e_el_out;

[a,b] = size(s_el_final);
l_period = b-mod(b,25);

s_el_final = s_el_final+adjust_el;
n_el_final = n_el_final+adjust_el;
% w_el_final = w_el_final-adjust_el;
e_el_final = e_el_final+adjust_el;

n_u_final = n_u_out;
s_u_final = s_u_out;
% w_u_final = w_u_out;
 e_u_final = e_u_out;

n_v_final = n_v_out;
s_v_final = s_v_out;
% w_v_final = w_v_out;
e_v_final = e_v_out;

n_u2d_final = n_u2d_out+tide_n_u;
s_u2d_final = s_u2d_out+tide_s_u;
%w_u2d_final = w_u2d_out+tide_w_u;
e_u2d_final = e_u2d_out+tide_e_u;


n_v2d_final = n_v2d_out+tide_n_v;
s_v2d_final = s_v2d_out+tide_s_v;
%w_v2d_final = w_v2d_out+tide_w_v;
e_v2d_final = e_v2d_out+tide_e_v;

% n_u2d_final = n_u2d_out;
% s_u2d_final = s_u2d_out;
% %w_u2d_final = w_u2d_out;
% e_u2d_final = e_u2d_out;
% 
% 
% n_v2d_final = n_v2d_out;
% s_v2d_final = s_v2d_out;
% %w_v2d_final = w_v2d_out;
% e_v2d_final = e_v2d_out;

n_temp = n_temp_out;
n_s = n_s_out;
s_temp = s_temp_out;
s_s = s_s_out;
e_temp = e_temp_out;
e_s = e_s_out;

 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  n_el_final = tide_n_el;
% s_el_final =  tide_s_el;
% % w_el_final = tide_w_el;
%  e_el_final =  tide_e_el;
%  n_v2d_final = tide_n_v;
% s_v2d_final = tide_s_v;
% %w_v2d_final = tide_w_v;
% e_v2d_final = tide_e_v;
% n_u2d_final = tide_n_u;
% s_u2d_final = tide_s_u;
% %w_u2d_final =tide_w_u;
% e_u2d_final = tide_e_u;
% 
% n_u_final = n_u_out*0;
% s_u_final = s_u_out*0;
% % w_u_final = w_u_out*0;
%  e_u_final = e_u_out*0;
% n_v_final = n_v_out*0;
% s_v_final = s_v_out*0;
% % w_v_final = w_v_out*0;
%  e_v_final = e_v_out*0;
 %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%n_temp = n_temp2;
%s_temp = s_temp2;
%e_temp = e_temp2;
%n_s = n_s2;
%s_s = s_s2;
%e_s = e_s2;


if(NCS>0)
    %mud
    for i=1:NCS
        eval(strcat('n_mud',sprintf('%02d',i),' = zeros(r,gn.N,length(mud_time));'));
        eval(strcat('s_mud',sprintf('%02d',i),' = zeros(r,gn.N,length(mud_time));'));
%         eval(strcat('w_mud',sprintf('%02d',i),' = zeros(c,gn.N,length(mud_time));'));
        eval(strcat('e_mud',sprintf('%02d',i),' = zeros(c,gn.N,length(mud_time));'));
    end
end

if(NNS>0)
    %sand
    for i=1:NNS
        eval(strcat('n_sand',sprintf('%02d',i),' = zeros(r,gn.N,length(mud_time));'));
        eval(strcat('s_sand',sprintf('%02d',i),' = zeros(r,gn.N,length(sand_time));'));
%         eval(strcat('w_sand',sprintf('%02d',i),' = zeros(c,gn.N,length(sand_time));'));
        eval(strcat('e_sand',sprintf('%02d',i),' = zeros(c,gn.N,length(sand_time));'));
    end
end

% write variable
ncwrite(fn,'zeta_time',zeta_time);
ncwrite(fn,'v2d_time',v2d_time);
if(NCS>0)
    ncwrite(fn,'mud_time',mud_time);
end
if(NNS>0)
    ncwrite(fn,'sand_time',sand_time);
end
ncwrite(fn,'v3d_time',v3d_time);
ncwrite(fn,'temp_time',temp_time);
ncwrite(fn,'salt_time',salt_time);

if(dye_flag==1)
    if(dye_flag==1)
        ncwrite(fn,'dye_time',dye_time);
    end
end

ncwrite(fn,'zeta_north',n_el_final);
ncwrite(fn,'zeta_south',s_el_final);
% ncwrite(fn,'zeta_west',w_el_final);
 ncwrite(fn,'zeta_east',e_el_final);

ncwrite(fn,'ubar_north',n_u2d_final);
ncwrite(fn,'ubar_south',s_u2d_final);
% ncwrite(fn,'ubar_west',w_u2d_final);
 ncwrite(fn,'ubar_east',e_u2d_final);

ncwrite(fn,'vbar_north',n_v2d_final);
ncwrite(fn,'vbar_south',s_v2d_final);
% ncwrite(fn,'vbar_west',w_v2d_final);
 ncwrite(fn,'vbar_east',e_v2d_final);

ncwrite(fn,'u_south',s_u_final);
ncwrite(fn,'u_north',n_u_final);
% ncwrite(fn,'u_west',w_u_final);
 ncwrite(fn,'u_east',e_u_final);

ncwrite(fn,'v_south',s_v_final);
ncwrite(fn,'v_north',n_v_final);
% ncwrite(fn,'v_west',w_v_final);
 ncwrite(fn,'v_east',e_v_final);

ncwrite(fn,'temp_north',n_temp);
ncwrite(fn,'temp_south',s_temp);
% ncwrite(fn,'temp_west',w_temp);
 ncwrite(fn,'temp_east',e_temp);

ncwrite(fn,'salt_north',n_s);
ncwrite(fn,'salt_south',s_s);
% ncwrite(fn,'salt_west',w_s);
 ncwrite(fn,'salt_east',e_s);

if(dye_flag==1)
    if(year==2021)
        n_tracer01 = zeros(r,gn.N,length(dye_time));
        s_tracer01 = zeros(r,gn.N,length(dye_time));
        e_tracer01 = zeros(c,gn.N,length(dye_time));
    end
 ncwrite(fn,'dye_north_01',n_tracer01);
 ncwrite(fn,'dye_south_01',s_tracer01);
 % ncwrite(fn,'dye_west_01',w_tracer01);
 ncwrite(fn,'dye_east_01',e_tracer01);
end

if(NCS>0)
    for i=1:NCS
    eval(strcat('ncwrite(fn,''mud_north_',sprintf('%02d',i),''',n_mud',sprintf('%02d',i),');'));
    eval(strcat('ncwrite(fn,''mud_south_',sprintf('%02d',i),''',s_mud',sprintf('%02d',i),');'));
%     eval(strcat('ncwrite(fn,''mud_west_',sprintf('%02d',i),''',w_mud',sprintf('%02d',i),');'));
    eval(strcat('ncwrite(fn,''mud_east_',sprintf('%02d',i),''',e_mud',sprintf('%02d',i),');'));
    end
end

if(NNS>0)
    for i=1:NNS
    eval(strcat('ncwrite(fn,''sand_north_',sprintf('%02d',i),''',n_sand',sprintf('%02d',i),');'));
    eval(strcat('ncwrite(fn,''sand_south_',sprintf('%02d',i),''',s_sand',sprintf('%02d',i),');'));
%     eval(strcat('ncwrite(fn,''sand_west_',sprintf('%02d',i),''',w_sand',sprintf('%02d',i),');'));
    eval(strcat('ncwrite(fn,''sand_east_',sprintf('%02d',i),''',e_sand',sprintf('%02d',i),');'));
    end
end


for i=1:l_period-24
    z_north(i) = mean(squeeze(n_el_final(447,i:i+24)));
    z_south(i) = mean(squeeze(s_el_final(200,i:i+24)));
%     z_west(i) = mean(squeeze(w_el_final(5,i:i+24)));
    z_east(i) = mean(squeeze(e_el_final(5,i:i+24)));
end

figure(1);
plot(z_north,'r');
hold on;
plot(z_south,'g');
hold on;
% plot(z_west,'y');
% hold on;
plot(z_east,'b');
hold on;

legend('North','South','East');








