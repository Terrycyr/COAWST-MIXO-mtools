clear all;
init_file = 'biotest_ini.nc'; delete(init_file);

grd = './bio1dtest_grd.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = 21;
N=gn.N;
Nbed = 0;
NNS = 0;
NCS =0;
Nveg=0;
NBT=0;
NDYE=1;
sed_flag=0;
dye_flag=0;
tidal_hycom = 1;
nonzero_ini = 0.;
el_adjust = 0.0;

create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,Nveg,NBT,NDYE,0);
const_initialize(init_file,0.)


[r,c] = size(lon);

Vtransform = ncread(grd,'Vtransform');                         
Vstretching = ncread(grd,'Vstretching');
theta_s = ncread(grd,'theta_s');                     
theta_b = ncread(grd,'theta_b');                     
Tcline = ncread(grd,'Tcline');
hc = ncread(grd,'hc');
Cs_w = ncread(grd,'Cs_w');
Cs_r = ncread(grd,'Cs_r');
s_w = ncread(grd,'s_w');
s_rho = ncread(grd,'s_rho');

settling_vel(1:r,1:c) = 1.0e-04;
erosion_stress(1:r,1:c) = 5.0e-05;
grain_diameter(1:r,1:c) = 8.0e-06;

spherical = ncread(grd,'spherical');
if(strcmp(spherical,'T'))
    ncwrite(init_file,'spherical',1);
else
    ncwrite(init_file,'spherical',0);
end

ncwrite(init_file,'theta_s',theta_s);
ncwrite(init_file,'theta_b',theta_b);
ncwrite(init_file,'Tcline',Tcline);
ncwrite(init_file,'Cs_r',Cs_r);
ncwrite(init_file,'Cs_w',Cs_w);
ncwrite(init_file,'s_w',s_w);
ncwrite(init_file,'s_rho',s_rho);
ncwrite(init_file,'hc',hc);
ncwrite(init_file,'Vtransform',Vtransform);
ncwrite(init_file,'Vstretching',Vstretching);
ncwrite(init_file,'spherical',spherical);


ncwrite(init_file,'grain_diameter',grain_diameter);
ncwrite(init_file,'settling_vel',settling_vel);
ncwrite(init_file,'erosion_stress',erosion_stress);

ncwrite(init_file,'ocean_time',0);

tmp = ncread(init_file,'temp');
tmp = ones(size(tmp)).*25;
ncwrite(init_file,'temp',tmp);
