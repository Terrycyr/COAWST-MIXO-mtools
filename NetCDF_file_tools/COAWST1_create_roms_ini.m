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
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\NC_file_generation');

year = 2002;
nontidal_dataset_dir = '../Non_tidal_component_preprocessing/HYCOM/';
nontidal_dataset_source = 'GLBv0.08/expt_53.X';

init_file = ['WFS_',num2str(year),'_ini.nc']; delete(init_file);
grd_name =  '../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N =length(ncread(grd_name,'Cs_r'));
N=gn.N;
Nbed = 0;
NNS = 0;
NCS = 0;
Nveg=0;
NBT=0;
NDYE=0;
sed_flag=0;
dye_flag=0;
tidal_hycom = 0;
nonzero_ini = 0;
el_adjust = 0.0;

mgl2kgm3 = 0.001;

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

settling_vel(1:r,1:c) = 9.8985e-04;
erosion_stress(1:r,1:c) = 6.8986e-05;
grain_diameter(1:r,1:c) = 3.2000e-05;

if(sed_flag==1)
    load(strcat('../Water_Atlas/','bay_ini_',num2str(year),'.mat'),'tss_bay');
    load('../Sediment/usSEABED/WFS_seabed.mat');

    %Parameters
    water_silt_clay_frac = 0.4;
    water_sand_frac = 0.6;

    if(NCS>0)
        mud_poros(1:r,1:c,1:Nbed) = 0.672d0;
    end
    if(NNS>0)
        sand_poros(1:r,1:c,1:Nbed) = 0.672d0;
    end

    grain_density(1:r,1:c) = 2650;
    bed_thickness(1:r,1:c,1) = 0.05;
    bed_thickness(1:r,1:c,2) = 0.05;
    bed_thickness(1:r,1:c,3) = 3.00;
    bed_age(1:r,1:c,1:Nbed) = 0.d0;

    if(NCS>0)
        bed_porosity = mud_poros;
    end
    if(NNS>0)
        bed_porosity = sand_poros;
    end
    bed_biodiff(1:r,1:c,1:Nbed) = 0.0;
    bed_tau_crit(1:r,1:c,1:Nbed) = 0.1;
    dmix_offset(1:r,1:c) = -0.4690;
    dmix_slope(1:r,1:c) = 1;
    dmix_time(1:r,1:c) = 0;
    ripple_height(1:r,1:c) = 0.01;
    
    %Initialize
    if(NCS>0)
        for i=1:NCS
            eval(strcat('mud_',sprintf('%02d',i),'(1:r,1:c,1:gn.N) = 0;'));
        end

        mud_01(1:r,1:c,1:gn.N) = tss_bay*mgl2kgm3*water_silt_clay_frac;
        mud_02(1:r,1:c,1:gn.N) = tss_bay*mgl2kgm3*water_sand_frac;

        %Remove NaNs
        for i=1:NCS
            eval(strcat('mud_',sprintf('%02d',i),'(isnan(mud_',sprintf('%02d',i),')) = 0;'));
        end

        for i=1:Nbed
            mudfrac_01(:,:,i) = silt_clay_frac;
            mudfrac_02(:,:,i) = sand_frac;

            %make sure none of classes is less than 5%
            mudfrac_02(mudfrac_01<0.05) = 0.95;
            mudfrac_01(mudfrac_01<0.05) = 0.05;

            mudfrac_01(mudfrac_02<0.05) = 0.95;
            mudfrac_02(mudfrac_02<0.05) = 0.05;


            for j=1:NCS
                eval(strcat('mudmass_',sprintf('%02d',j),'(:,:,i)',...
                    ' = bed_thickness(:,:,i).*grain_density.*(1-mud_poros(:,:,i)).*mudfrac_',...
                    sprintf('%02d',j),'(:,:,i);'));

                eval(strcat('mudfrac_',sprintf('%02d',j),'(isnan(mudfrac_',sprintf('%02d',j),')) = 0.0;'));
                eval(strcat('mudmass_',sprintf('%02d',j),'(isnan(mudmass_',sprintf('%02d',j),')) = 0.0;'));
            end
        end
    end

    if(NNS>0)
        for i=1:NNS
            eval(strcat('sand_',sprintf('%02d',i),'(1:r,1:c,1:gn.N) = 0;'));
        end

        sand_01(1:r,1:c,1:gn.N) = tss_bay*water_silt_clay_frac;
        sand_02(1:r,1:c,1:gn.N) = tss_bay*water_sand_frac;

        %Remove NaNs
        for i=1:NNS
            eval(strcat('sand_',sprintf('%02d',i),'(isnan(sand_',sprintf('%02d',i),')) = 0;'));
        end

        for i=1:Nbed
            sandfrac_01(:,:,i) = silt_clay_frac(:,:,i);
            sandfrac_02(:,:,i) = sand_frac(:,:,i);


            for j=1:NNS
                eval(strcat('sandmass_',sprintf('%02d',j),'(:,:,i)',...
                    ' = bed_thickness(:,:,i).*grain_density.*(1-sand_poros(:,:,i)).*sandfrac_',...
                    sprintf('%02d',j),'(:,:,i);'));

                eval(strcat('sandfrac_',sprintf('%02d',j),'(isnan(sandfrac_',sprintf('%02d',j),')) = 0.0;'));
                eval(strcat('sandmass_',sprintf('%02d',j),'(isnan(sandmass_',sprintf('%02d',j),')) = 0.0;'));
            end
        end
    end
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
    if(NCS>0)
        for i=1:NCS
            eval(strcat('ncwrite(init_file,''mud_',sprintf('%02d',i),''',mud_',sprintf('%02d',i),');'));
            eval(strcat('ncwrite(init_file,''mudfrac_',sprintf('%02d',i),''',mudfrac_',sprintf('%02d',i),');'));
            eval(strcat('ncwrite(init_file,''mudmass_',sprintf('%02d',i),''',mudmass_',sprintf('%02d',i),');'));
        end
    end

    if(NNS>0)
        for i=1:NNS
            eval(strcat('ncwrite(init_file,''sand_',sprintf('%02d',i),''',sand_',sprintf('%02d',i),');'));
            eval(strcat('ncwrite(init_file,''sandfrac_',sprintf('%02d',i),''',sandfrac_',sprintf('%02d',i),');'));
            eval(strcat('ncwrite(init_file,''sandmass_',sprintf('%02d',i),''',sandmass_',sprintf('%02d',i),');'));
        end
    end

    ncwrite(init_file,'bed_thickness',bed_thickness);
    ncwrite(init_file,'bed_age',bed_age);
    ncwrite(init_file,'bed_porosity',bed_porosity);
    ncwrite(init_file,'bed_biodiff',bed_biodiff);
    ncwrite(init_file,'bed_tau_crit',bed_tau_crit);
    
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

