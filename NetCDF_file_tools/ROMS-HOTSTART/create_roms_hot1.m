clear all;
% toolbox path
 addpath(path,'../COAWST_v3.4/Tools/mfiles/rutgers/utility');
 addpath(path,'../COAWST_v3.4/Tools/mfiles/roms_clm');

init_file = 'WFS_2011_hot.nc';
grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v2.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
angle = ncread(grd_name,'angle');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N = 20;
N=gn.N;
Nbed = 0;
NNS = 0;
NCS =0;
Nveg=0;
sed_flag=0;
create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,Nveg);

[r,c] = size(lon);

theta_b = 0;
theta_s = 0;
Tcline = 0;
hmin=0.5;
hc=min([hmin,Tcline]);

settling_vel(1:r,1:c) = 1.0e-04;
erosion_stress(1:r,1:c) = 5.0e-05;
grain_diameter(1:r,1:c) = 8.0e-06;

his = './WFS_2011_his_0365.nc';

s_hot = ncread(his,'salt',[1 1 1 24], [Inf Inf Inf 1]);
temp_hot = ncread(his,'temp',[1 1 1 24], [Inf Inf Inf 1]);
zeta_hot = ncread(his,'zeta',[1 1 24], [Inf Inf 1]);
u = ncread(his,'u_eastward',[1 1 1 24], [Inf Inf Inf 1]);
v = ncread(his,'v_northward',[1 1 1 24], [Inf Inf Inf 1]);
ubar = ncread(his,'ubar_eastward',[1 1 24], [Inf Inf 1]);
vbar = ncread(his,'vbar_northward',[1 1 24], [Inf Inf 1]);

s_hot(isnan(s_hot)) = 0;
temp_hot(isnan(temp_hot)) = 0;
zeta_hot(isnan(zeta_hot)) = 0;
u(isnan(u)) = 0;
v(isnan(v)) = 0;
ubar(isnan(ubar)) = 0;
vbar(isnan(vbar)) = 0;

u2d_hot = rho2u_2d_mw(cos(angle).*ubar+sin(angle).*vbar);
v2d_hot = rho2v_2d_mw(-sin(angle).*ubar+cos(angle).*vbar);


for i=1:N
    u_hot(:,:,i) = rho2u_2d_mw(cos(angle).*u(:,:,i)+sin(angle).*v(:,:,i));
    v_hot(:,:,i) = rho2v_2d_mw(-sin(angle).*u(:,:,i)+cos(angle).*v(:,:,i));
end
[r,c] = size(lon);

spherical = ncread(grd_name,'spherical');
if(strcmp(spherical,'T'))
    ncwrite(init_file,'spherical',1);
else
    ncwrite(init_file,'spherical',0);
end

ncwrite(init_file,'Vtransform',2);
ncwrite(init_file,'Vstretching',4);
ncwrite(init_file,'theta_b',theta_b);
ncwrite(init_file,'theta_s',theta_s);
ncwrite(init_file,'Tcline',Tcline);
ncwrite(init_file,'hc',hc);

if(sed_flag==1)
    ncwrite(init_file,'mud_01',mud_01);
    ncwrite(init_file,'mud_02',mud_02);
    ncwrite(init_file,'mud_03',mud_03);
    ncwrite(init_file,'mud_04',mud_04);
    ncwrite(init_file,'mud_05',mud_05);

    ncwrite(init_file,'mudfrac_01',mudfrac_01);
    ncwrite(init_file,'mudfrac_02',mudfrac_02);
    ncwrite(init_file,'mudfrac_03',mudfrac_03);
    ncwrite(init_file,'mudfrac_04',mudfrac_04);
    ncwrite(init_file,'mudfrac_05',mudfrac_05);

    ncwrite(init_file,'mudmass_01',mudmass_01);
    ncwrite(init_file,'mudmass_02',mudmass_02);
    ncwrite(init_file,'mudmass_03',mudmass_03);
    ncwrite(init_file,'mudmass_04',mudmass_04);
    ncwrite(init_file,'mudmass_05',mudmass_05);

    ncwrite(init_file,'bed_thickness',bed_thickness);
    ncwrite(init_file,'bed_age',bed_age);
    ncwrite(init_file,'bed_porosity',bed_porosity);
    ncwrite(init_file,'bed_biodiff',bed_biodiff);
    
    ncwrite(init_file,'ripple_height',ripple_height);
    ncwrite(init_file,'dmix_offset',dmix_offset);
    ncwrite(init_file,'dmix_slope',dmix_slope);
    ncwrite(init_file,'dmix_time',dmix_time);
    ncwrite(init_file,'grain_density',grain_density);
end

ncwrite(init_file,'grain_diameter',grain_diameter);
ncwrite(init_file,'settling_vel',settling_vel);
ncwrite(init_file,'erosion_stress',erosion_stress);

ncwrite(init_file,'ocean_time',0);
ncwrite(init_file,'salt',s_hot);
ncwrite(init_file,'temp',temp_hot);
ncwrite(init_file,'u',u_hot);
ncwrite(init_file,'ubar',u2d_hot);
ncwrite(init_file,'v',v_hot);
ncwrite(init_file,'vbar',v2d_hot);
ncwrite(init_file,'zeta',zeta_hot);


