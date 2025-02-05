function [kb_ini,kb_lon,kb_lat,kb_dat] = gen_kb_ini(t_min,t_max,kb_bg,mat_data,method)

    %Toolbox
    addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
    addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');

    %Params.
    grd = 'C:\Users\cheny\Desktop\EcoHAB\Model_grid\ROMS_WFS_new.nc';
    bay_mask = ncread('C:\Users\cheny\Desktop\EcoHAB\Model_grid\ROMS_WFS_new_bay_mask.nc','mask_rho');

    lon = ncread(grd,'lon_rho');
    lat = ncread(grd,'lat_rho');
    mask = ncread(grd,'mask_rho');


    [kb_lon,kb_lat, kb_dep, kb_dat, kb_dep_all, kb_dat_all] = get_FWC_kb(t_min,t_max,mat_data,method);

    filter_pos = find(kb_lat>29);
    kb_lon(filter_pos) = [];
    kb_lat(filter_pos) = [];
    kb_dat(filter_pos) = [];
    
    kb_ini = griddata(kb_lon(kb_dat>kb_bg),kb_lat(kb_dat>kb_bg),kb_dat(kb_dat>kb_bg),lon,lat,'natural');
    kb_ini(mask==0) = 0;
    kb_ini(isnan(kb_ini)) = 0;

%     figure;
%     m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
%     m_contourf(lon,lat,kb_ini,100,'linestyle','none');
%     hold on;
%     m_gshhs_h('patch',[.6 .6 .6]);
%     hold on;
%     m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
%     colorbar;
%     set(gcf,'color', [1 1 1]);
%     xlabel('Longitude','fontsize',9);
%     ylabel('Latitude','fontsize',9);

end
