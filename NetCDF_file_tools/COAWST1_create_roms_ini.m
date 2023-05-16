clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check following variables when setting up the script
% 1. Grid name.
% 2. Year
% 3. nontidal_dataset_dir and nontidal_dataset_source
% 4. Output file name
% 5. sed_flag, dye_flag and nonzero_ini
% 6. tidal_hycom and el_adjust, mostly tidal_hycom = 0 to use non-tidal
% parts from HYCOM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

year = 2021;
nontidal_dataset_dir = '../Non_tidal_component_preprocessing/HYCOM/Detided_GOM_HYCOM/';
nontidal_dataset_source = 'GOMu0.04/expt_90.1m000';
%nontidal_dataset_source = 'GLBy0.08/expt_93.0';

init_file = ['WFS_',num2str(year),'_ini_detide.nc']; delete(init_file);
%grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
grd_name =  '../Model_grid/ROMS_WFS_Piney.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N =length(ncread(grd_name,'Cs_r'));
N=gn.N;
Nbed = 0;
NNS = 0;
NCS =0;
Nveg=0;
NBT=0;
NDYE=1;
sed_flag=0;
dye_flag=1;
tidal_hycom = 1;
nonzero_ini = 0;
el_adjust = 0.0;

load(strcat('../Tidal_component_preprocessing/tide_ini_',num2str(year),'.mat'));
load(strcat(nontidal_dataset_dir,'non_tidal_ini_',num2str(year),'.mat'));

if(dye_flag==1)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     load(strcat('../GOM_preprocessing/ocean_nutrients_ini_',num2str(year),'.mat'));
    load(strcat('../Water_Atlas/','bay_ini_',num2str(year),'.mat'));


%     for i=1:gn.N
%         tmp = no3_roms(:,:,i);
%         tmp(~isnan(no23_bay(:,:,i))) = no23_bay(~isnan(no23_bay(:,:,1)));
%         NO23(:,:,i) = tmp;
%     end
% 
%     tracer01 = NO23;
  
    tracer01 = zeros(size(sal_bay));
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else
    load(strcat('../Water_Atlas/','bay_ini_',num2str(year),'.mat'),'sal_bay');
end

create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,Nveg,NBT,NDYE,0, nontidal_dataset_source);

[r,c] = size(lon);
N= length(ncread(grd_name,'Cs_r'));

Vtransform = ncread(grd_name,'Vtransform');                         
Vstretching = ncread(grd_name,'Vstretching');
theta_s = ncread(grd_name,'theta_s');                     
theta_b = ncread(grd_name,'theta_b');                     
Tcline = ncread(grd_name,'Tcline');
hc = ncread(grd_name,'hc');
Cs_w = ncread(grd_name,'Cs_w');
Cs_r = ncread(grd_name,'Cs_r');
s_w = ncread(grd_name,'s_w');
s_rho = ncread(grd_name,'s_rho');

settling_vel(1:r,1:c) = 1.0e-04;
erosion_stress(1:r,1:c) = 5.0e-05;
grain_diameter(1:r,1:c) = 8.0e-06;

if(sed_flag==1)
    load('./clay_grid2_zh.mat');
    load('./silt_grid2_zh.mat');
    load('./sand_grid2_zh.mat');
    mud_01(1:r,1:c,gn.N) = 0;
    mud_02(1:r,1:c,gn.N) = 0;
    mud_03(1:r,1:c,gn.N) = 0;
    mud_04(1:r,1:c,gn.N) = 0;
    mud_05(1:r,1:c,gn.N) = 0;

    mud_poros(1:r,1:c,1:Nbed) = 0.672d0;
    bed_thickness(1:r,1:c,1) = 0.05; 
    bed_thickness(1:r,1:c,2) = 0.03; 
    bed_thickness(1:r,1:c,3) = 3.00; 

    grain_density(1:r,1:c) = 2650;
    
    for i=1:Nbed
        mudfrac_01(:,:,i) = clay_grid2;
        mudfrac_02(:,:,i) = silt_grid2;
        mudfrac_03(:,:,i) = sand_grid2;
        mudfrac_04(:,:,i) = zeros(size(clay_grid2));
        mudfrac_05(:,:,i) = zeros(size(clay_grid2));

        mudmass_01(:,:,i) = bed_thickness(:,:,i).*grain_density.*(1-mud_poros(:,:,i)).*mudfrac_01(:,:,i);
        mudmass_02(:,:,i) = bed_thickness(:,:,i).*grain_density.*(1-mud_poros(:,:,i)).*mudfrac_02(:,:,i);
        mudmass_03(:,:,i) = bed_thickness(:,:,i).*grain_density.*(1-mud_poros(:,:,i)).*mudfrac_03(:,:,i);
        mudmass_04(:,:,i) = bed_thickness(:,:,i).*grain_density.*(1-mud_poros(:,:,i)).*mudfrac_04(:,:,i);
        mudmass_05(:,:,i) = bed_thickness(:,:,i).*grain_density.*(1-mud_poros(:,:,i)).*mudfrac_05(:,:,i);

        mudfrac_01(isnan(mudfrac_01)) = 0.0;
        mudfrac_02(isnan(mudfrac_02)) = 0.0;
        mudfrac_03(isnan(mudfrac_03)) = 0.0;
        mudfrac_04(isnan(mudfrac_03)) = 0.0;
        mudfrac_05(isnan(mudfrac_03)) = 0.0;

        mudmass_01(isnan(mudmass_01)) = 0.0;
        mudmass_02(isnan(mudmass_02)) = 0.0;
        mudmass_03(isnan(mudmass_03)) = 0.0;
        mudmass_04(isnan(mudmass_03)) = 0.0;
        mudmass_05(isnan(mudmass_03)) = 0.0;
    end

    bed_age(1:r,1:c,1:Nbed) = 0.d0;
    bed_porosity = mud_poros;
    bed_biodiff(1:r,1:c,1:Nbed) = 0.0;
    bed_tau_crit(1:r,1:c,1:Nbed) = 0.05;
    dmix_offset(1:r,1:c) = -0.4690;
    dmix_slope(1:r,1:c) = 1;
    dmix_time(1:r,1:c) = 0;
    ripple_height(1:r,1:c) = 0.01;
end    

if(tidal_hycom==1)
    el_final = el_ini_out+el_adjust;
    u2d_final = u2d_ini_out;
    v2d_final = v2d_ini_out;
else
    el_final = el_ini_out+tide_el_ini+el_adjust;
    u2d_final = u2d_ini_out+tide_u_ini;
    v2d_final = v2d_ini_out+tide_v_ini;
end

u_final = u_ini_out;
v_final = v_ini_out;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
el_final = el_final*nonzero_ini;
u2d_final = u2d_final*nonzero_ini;
v2d_final = v2d_final*nonzero_ini;
u_final = u_final*nonzero_ini;
v_final = v_final*nonzero_ini;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

spherical = ncread(grd_name,'spherical');
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

% el_final = el_final;
% u_final = u_final;
% v_final = v_final;
% u2d_final = u2d_final;
% v2d_final =v2d_final;
%s_ini_out = ones(r,c,N)*36;
%temp_ini_out = ones(r,c,N)*24;
%temp_ini_out = temp_roms;
%s_ini_out  = sal_roms;

%temp_ini_out = imgaussfilt(temp_ini_out,1.2);
%s_ini_out = imgaussfilt(s_ini_out,1.2);

%for i=1:N
%    tmp = s_ini_out(:,:,i);
%    tmp(dep>500&tmp>34.7) = 36;
%end

k=0;
for i=N:-1:1
    k=k+1;
    tmp = temp_ini_out(:,:,i);
    tmp(mask==0) =0;
    temp_ini_out(:,:,i) = tmp;
    
    tmp = s_ini_out(:,:,i);
    tmp(mask==0) =0;
    s_ini_out(:,:,i) = tmp;  

    tmp = s_ini_out(:,:,i);
    tmp2 = sal_bay(:,:,k);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    s_ini_out(:,:,i) = tmp;
end

[ru,cu] = size(u_final);
[rv,cv] = size(v_final);
ncwrite(init_file,'ocean_time',0);
ncwrite(init_file,'salt',s_ini_out);
ncwrite(init_file,'temp',temp_ini_out);
ncwrite(init_file,'u',u_final);
ncwrite(init_file,'ubar',u2d_final);
ncwrite(init_file,'v',v_final);
ncwrite(init_file,'vbar',v2d_final);
ncwrite(init_file,'zeta',el_final);

if(dye_flag==1)
    ncwrite(init_file,'dye_01',tracer01);
end

