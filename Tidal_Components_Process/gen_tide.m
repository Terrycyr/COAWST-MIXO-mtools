% By Yuren Chen, 2020-07-05, Version 1.0
%--------------------------------------------------------------------------
% Program generate tidal boudaries from  output files of 
% OSU Tidal Prediction Software.
% Due to the mismatch of land masks and bathymetry between OTPS
% dataset and ROMS model grid, extrapolation and flux correction are
% performed here. 
%--------------------------------------------------------------------------

clear all; close all;

% toolbox path
addpath(path,'../../COAWST/Tools/mfiles/roms_clm');

% model grid
fn = '../Model_grid/ROMS_WFS_new.nc';

% Output data length 
dt=1/24;
d1=datenum(2003,1,1,0,0,0);
d2=datenum(2003,12,31,24,0,0);
tide_out_date = d1:dt:d2;
year = datevec(d1);
year = year(1);
%Boundary used
bnd_flag = [1 1 0 1]; %E S W N

% OTPS output
fn1=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_el_e.out'];
fn2=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_el_s.out'];
fn3=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_el_w.out'];
fn4=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_el_n.out'];

fn5=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_u_e.out'];
fn6=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_u_s.out'];
fn7=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_u_w.out'];
fn8=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_u_n.out'];

fn9=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_v_e.out'];
fn10=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_v_s.out'];
fn11=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_v_w.out'];
fn12=['./OTPS/COAWST_TIDE_OUTPUT_',num2str(year),'/coawst_v_n.out'];

%--------------------- Read from model grid--------------------------------
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

save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_out_date','-v7.3');
%EAST
if(bnd_flag(1)==1)
    grd.blon = lon(end,:)';
    grd.blat = lat(end,:)';
    grd.bmask = mask(end,:)';
    grd.bmask_u = mask_u(end,:)';
    grd.bmask_v = mask_v(end,:)';
    grd.bangle = angle(end,:)';
    grd.bangle_u = angle_u(end,:)';
    grd.bangle_v = angle_v(end,:)';
    grd.blon_u = lon_u(end,:)';
    grd.blat_u = lat_u(end,:)';
    grd.blon_v = lon_v(end,:)';
    grd.blat_v = lat_v(end,:)';
    grd.bh = h(end,:)';
    grd.bh_u = h_u(end,:)';
    grd.bh_v = h_v(end,:)';
    
    e_el = get_OTPS_el_result(grd, tide_out_date,fn1);
    [e_u, e_v, e_u4, e_v4] = get_OTPS_uv_result(grd, tide_out_date,fn5,fn9);
    
    % Figures to help checking
    figure;
    quiver(grd.blon_v,grd.blat_v,e_u(1:end-1,1,1),e_v(:,1,1));
    hold on;
    quiver(grd.blon_v,grd.blat_v,e_u4(1:end-1,1,1),e_v4(:,1,1));
    legend('After','Before');

    tide_e_el = e_el;
    tide_e_u = e_u;
    tide_e_v = e_v;
    
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_e_el','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_e_u','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_e_v','-append');
end

%SOUTH
if(bnd_flag(2)==1)
    grd.blon = lon(:,1);
    grd.blat = lat(:,1);
    grd.bmask = mask(:,1);
    grd.bmask_u = mask_u(:,1);
    grd.bmask_v = mask_v(:,1);
    grd.bangle = angle(:,1);
    grd.bangle_u = angle_u(:,1);
    grd.bangle_v = angle_v(:,1);
    grd.blon_u = lon_u(:,1);
    grd.blat_u = lat_u(:,1);
    grd.blon_v = lon_v(:,1);
    grd.blat_v = lat_v(:,1);
    grd.bh = h(:,1);
    grd.bh_u = h_u(:,1);
    grd.bh_v = h_v(:,1);
    
    s_el = get_OTPS_el_result(grd, tide_out_date,fn2);
    [s_u, s_v, s_u4, s_v4] = get_OTPS_uv_result(grd, tide_out_date,fn6,fn10);
    
    figure;
    quiver(grd.blon_u,grd.blat_u,s_u(:,1,1),s_v(1:end-1,1,1));
    hold on;
    quiver(grd.blon_u,grd.blat_u,s_u4(:,1,1),s_v4(1:end-1,1,1));
    legend('After','Before');
    
    tide_s_el = s_el;
    tide_s_u = s_u;
    tide_s_v = s_v;

    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_s_el','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_s_u','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_s_v','-append');
end

%WEST
if(bnd_flag(3)==1)
    grd.blon = lon(1,:)';
    grd.blat = lat(1,:)';
    grd.bmask = mask(1,:)';
    grd.bmask_u = mask_u(1,:)';
    grd.bmask_v = mask_v(1,:)';
    grd.bangle = angle(1,:)';
    grd.bangle_u = angle_u(1,:)';
    grd.bangle_v = angle_v(1,:)';
    grd.blon_u = lon_u(1,:)';
    grd.blat_u = lat_u(1,:)';
    grd.blon_v = lon_v(1,:)';
    grd.blat_v = lat_v(1,:)';
    grd.bh = h(1,:)';
    grd.bh_u = h_u(1,:)';
    grd.bh_v = h_v(1,:)';
    
    w_el = get_OTPS_el_result(grd, tide_out_date,fn3);
    [w_u, w_v, w_u4, w_v4] = get_OTPS_uv_result(grd, tide_out_date,fn7,fn11);
    
    figure;
    quiver(grd.blon_v',grd.blat_v',w_u(1:end-1,1,1),w_v(:,1,1));
    hold on;
    quiver(grd.blon_v',grd.blat_v',w_u4(1:end-1,1,1),w_v4(:,1,1));
    legend('After','Before');
    
    tide_w_el = w_el;
    tide_w_u = w_u;
    tide_w_v = w_v;

    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_w_el','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_w_u','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_w_v','-append');
end

%NORTH
if(bnd_flag(4)==1)
    grd.blon = lon(:,end);
    grd.blat = lat(:,end);
    grd.bmask = mask(:,end);
    grd.bmask_u = mask_u(:,end);
    grd.bmask_v = mask_v(:,end);
    grd.bangle = angle(:,end);
    grd.bangle_u = angle_u(:,end);
    grd.bangle_v = angle_v(:,end);
    grd.blon_u = lon_u(:,end);
    grd.blat_u = lat_u(:,end);
    grd.blon_v = lon_v(:,end);
    grd.blat_v = lat_v(:,end);
    grd.bh = h(:,end);
    grd.bh_u = h_u(:,end);
    grd.bh_v = h_v(:,end);
    
    n_el = get_OTPS_el_result(grd, tide_out_date,fn4);
    [n_u, n_v, n_u4, n_v4] = get_OTPS_uv_result(grd, tide_out_date,fn8,fn12);
    
    figure;
    quiver(grd.blon_u,grd.blat_u,n_u(:,1,1),n_v(1:end-1,1,1));
    hold on;
    quiver(grd.blon_u,grd.blat_u,n_u4(:,1,1),n_v4(1:end-1,1,1));
    legend('After','Before');

    tide_n_el = n_el;
    tide_n_u = n_u;
    tide_n_v = n_v;
    
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_n_el','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_n_u','-append');
    save(strcat('tide_bnd_',num2str(year),'.mat'),'tide_n_v','-append');
end


    