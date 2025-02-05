clear all; close all;

grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');

N= length(ncread(grd_name,'Cs_r'));
n_layer = N;
Vtransform = ncread(grd_name,'Vtransform');                         
Vstretching = ncread(grd_name,'Vstretching');
THETA_S = ncread(grd_name,'theta_s');                     
THETA_B = ncread(grd_name,'theta_b');                     
TCLINE = ncread(grd_name,'Tcline');
hc = ncread(grd_name,'hc');

% Calculate Sigma for each layer
sc = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,1,0,0);
sw = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,5,1,0,0);


ncwrite(grd_name,'Cs_r',squeeze(sc));
ncwrite(grd_name,'Cs_w',squeeze(sw));


