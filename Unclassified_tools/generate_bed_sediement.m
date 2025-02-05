clear all; close all;

%Toolbox
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');

%Params.
grd = '../Model_grid/ROMS_WFS_Piney.nc';
mask = ncread(grd,'mask_rho');
bay_mask = ncread('../Model_grid/ROMS_WFS_Piney_bay_mask.nc','mask_rho');
year = 2022;
Nbed = 3;

lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
h = ncread(grd,'h');

[r,c] =size(lon);

%Assume coarser sediment in areas further to Tampa Bay and Charlotte Harbour. 
scale_distance1 = get_scale_distance(lon,lat,bay_mask,mask,0.1,0.001);

%Constant initials
silt_clay_frac0 = 0.57*ones(r,c,Nbed);

for i  = 1:Nbed
    silt_clay_frac(:,:,i) = silt_clay_frac0(:,:,i).*scale_distance1;
    sand_frac(:,:,i) = 1-silt_clay_frac(:,:,i);
end

figure('units','pixels','position',[500 150 450 410]);
m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
hp = m_contourf(lon,lat,silt_clay_frac(:,:,1),100,'linestyle','none');
hold on;
m_gshhs_h('patch',[.6 .6 .6]);
hold on;
m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
colorbar;
set(gcf,'color', [1 1 1]);
xlabel('Longitude','fontsize',9);
ylabel('Latitude','fontsize',9); 

save(strcat('bed_sediment_ini_',num2str(year),'.mat'),'silt_clay_frac','sand_frac');



