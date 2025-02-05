% By Yuren Chen, 2021-04-07, Version 2.0
%--------------------------------------------------------------------------
% Program generate tidal initial field from 
% output files of OSU Tidal Prediction Software.
% Due to the mismatch of land masks and bathymetry between OTPS
% dataset and ROMS model grid, extrapolation and flux correction are
% performed here. 
%--------------------------------------------------------------------------

clear all; close all;

year = 2003;

% toolbox path
addpath(path,'../../COAWST/Tools/mfiles/roms_clm');

% model grid
%fn = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
fn = '../Model_grid/ROMS_WFS_new.nc';

% OTPS output
fn1=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_el_ini.out'];
fn4=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_u_ini.out'];
fn7=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_v_ini.out'];

fid1 = fopen(fn1,'r');
fid4 = fopen(fn4,'r');
fid7 = fopen(fn7,'r');

lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
mask = ncread(fn,'mask_rho');
mask_u = ncread(fn,'mask_u');
mask_v = ncread(fn,'mask_v');
lon_u = ncread(fn,'lon_u');
lon_v = ncread(fn,'lon_v');
lat_u = ncread(fn,'lat_u');
lat_v = ncread(fn,'lat_v');
h = ncread(fn,'h');
h_u = rho2u_2d_mw(h);
h_v = rho2v_2d_mw(h);
angle = ncread(fn,'angle');
[r,c] = size(mask);
[r1,c1] = size(lon_u);
[r2,c2] = size(lon_v);

for i= 1:r1
    angle_u(i,:) = (angle(i,:)+angle(i+1,:))/2;
end

for i= 1:c2
    angle_v(:,i) = (angle(:,i)+angle(:,i+1))/2;
end
   

nday = 0;

for i=1:6
    dat1=fgetl(fid1);
    dat4=fgetl(fid4);
    dat7=fgetl(fid7);
end


for i=1:r
    for j=1:c
        dat1=fgetl(fid1);
        if(strcmp(dat1(22),'*'))
            el_ini(i,j) = 0;
            dep(i,j) = h(i,j);
        else
            el_ini(i,j) = str2num(dat1(46:56));
            dep(i,j) = str2num(dat1(57:end));
        end
    end
end

for i=1:r1
    for j=1:c1
        dat4=fgetl(fid4);
        if(strcmp(dat4(22),'*'))
            u_ini(i,j) = 0;
            uv_ini(i,j) = 0;
            dep_u(i,j) = h_u(i,j);
        else
            u_ini(i,j) = str2num(dat4(66:76))/100;
            uv_ini(i,j) = str2num(dat4(77:86))/100;
            dep_u(i,j) = str2num(dat4(87:end));
        end 
    end
end

for i=1:r2
    for j=1:c2
        dat7=fgetl(fid7);
        if(strcmp(dat7(22),'*'))
            v_ini(i,j) = 0;
            vu_ini(i,j) = 0;
            dep_v(i,j) = h_v(i,j);
        else
            v_ini(i,j) = str2num(dat7(77:86))/100;
            vu_ini(i,j) = str2num(dat7(66:76))/100;
            dep_v(i,j) = str2num(dat7(87:end));
        end 
    end
end

u4_ini = u_ini;
v4_ini = v_ini;

%rotate u,v to x,e
u_ini = cos(angle_u).*u4_ini+sin(angle_u).*uv_ini;
v_ini = -sin(angle_v).*vu_ini+cos(angle_v).*v4_ini;

%flux correction
cor_u = dep_u./h_u;
cor_v = dep_v./h_v;

cor_u(h_u<20) = 1;
cor_v(h_v<20) = 1;

u_ini = u_ini.*cor_u;
v_ini = v_ini.*cor_v;

el_ini(mask==0) = 0;
u_ini(mask_u==0) = 0;
v_ini(mask_v==0) = 0;

tide_el_ini = el_ini;
tide_u_ini = u_ini;
tide_v_ini = v_ini;

save(strcat('tide_ini_',num2str(year),'.mat'),'tide_el_ini','-v7.3');
save(strcat('tide_ini_',num2str(year),'.mat'),'tide_u_ini','-append');
save(strcat('tide_ini_',num2str(year),'.mat'),'tide_v_ini','-append');

