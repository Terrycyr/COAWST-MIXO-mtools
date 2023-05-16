clear all; close all;
addpath(path,genpath('./rutgers'));
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');

gn_coarse = '../../Model_grid/ROMS_WFS_Piney.nc';
gn_fine = './ROMS_TBays_100m.nc';
gn_contact = 'ROMS_TBays_contact.nc';
%All bays
%irange = [206 392];

delete(gn_fine);
delete(gn_contact);

%Tampa
irange = [206 288];
jrange = [85 158];
refine_factor = 3;

plot_interval = 1;

lon = ncread(gn_coarse,'lon_rho');
lat = ncread(gn_coarse,'lat_rho');
lon_psi = ncread(gn_coarse,'lon_psi');
lat_psi = ncread(gn_coarse,'lat_psi');
mask_psi = ncread(gn_coarse,'mask_psi');
dep = ncread(gn_coarse,'h');
mask = ncread(gn_coarse,'mask_rho');
bay_mask = ncread('../../Model_grid/ROMS_WFS_Piney_bay_mask.nc','mask_rho');


if(refine_factor==1)
    S = grid_extract(gn_coarse,gn_fine,irange(1),irange(2),jrange(1),jrange(2));
else
    S = coarse2fine(gn_coarse,gn_fine,refine_factor,irange(1),irange(2),jrange(1),jrange(2));
end


tmp = GridBuilder;
disp('Press any key after finnishing editing the mask: ')
pause;


[S2, G2] = contact({gn_coarse,gn_fine},gn_contact);

figure('units','pixels','position',[500 150 500 500]);
lon_min = min(min(S.lon_psi));
lon_max = max(max(S.lon_psi));
lat_min = min(min(S.lat_psi));
lat_max = max(max(S.lat_psi));
x0 = lon(1:plot_interval:end,1:plot_interval:end);
y0 = lat(1:plot_interval:end,1:plot_interval:end);
mask0 = mask(1:plot_interval:end,1:plot_interval:end);
x1 = S.lon_rho(1:plot_interval:end,1:plot_interval:end);
y1 = S.lat_rho(1:plot_interval:end,1:plot_interval:end);
m_proj('miller','lat',[lat_min lat_max],'long',[lon_min lon_max]);
hold on;
m_grid('linestyle','none','tickstyle','dd','fontsize',11);
hold on;
m_pcolor(x0,y0,mask0);
hold on;
m_mesh(x1,y1,ones(size(x1)),'EdgeColor','g','FaceColor','none');
hold on;
%m_gshhs_i('patch',[1 1 1],'edgecolor','b');