function [kb_ini,syn_ini,sal_ini,temp_ini,dpo4_ini,dnh4_ini,dnox_ini,don_ini,dop_ini,tn_ini,tp_ini,sio2_ini] = get_OLD_ECOHAB_ini(year)

    %Toolbox
    addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
    addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');

    %Params.
    grd = 'C:\Users\cheny\Desktop\EcoHAB\Model_grid\ROMS_WFS_New.nc';
    bay_mask = ncread('C:\Users\cheny\Desktop\EcoHAB\Model_grid\ROMS_WFS_Piney_bay_mask.nc','mask_rho');

    lon = ncread(grd,'lon_rho');
    lat = ncread(grd,'lat_rho');
    mask = ncread(grd,'mask_rho');

    [~,~,raw] = xlsread('C:\Users\cheny\Desktop\EcoHAB\Cruise_data\ECOHAB-Karenia Data -4 blooms from Access-WORKING COMBO FILE (ALL STATIONS) 3-31-23.xlsx',2);

    k=0;

    for i=4:764
        tmp = datenum(raw{i,8},'mm/dd/yyyy');
        tmp2 = datevec(tmp);
        if(sum(tmp2(1)==year)>0)
            k=k+1;
            dat_date(k) = datenum(raw{i,8},'mm/dd/yyyy');
            dat_lat(k) = raw{i,5};
            dat_lon(k) = raw{i,6};
            dat_dep(k) = raw{i,7};
            dat_kb(k) = raw{i,26}; %cells/L
            dat_syn(k) = raw{i,124}*1000; %cells/L
            dat_sal(k) = raw{i,43};
            dat_temp(k) = raw{i,45};
            
            dat_dpo4(k) = raw{i,64}*31/1000; %mgP/L
            dat_dnh4(k) = raw{i,65}*14/1000; %mgN/L
            dat_dnox(k) = raw{i,65}*14/1000; %mgN/L
            dat_don(k) = raw{i,82}*14/1000; %mgN/L
            dat_dop(k) = raw{i,83}*31/1000; %mgP/L
            dat_tn(k) = raw{i,84}*14/1000; %mgN/L
            dat_tp(k) = raw{i,85}*31/1000; %mgP/L
            dat_sio2(k) = raw{i,69}*28/1000; %mgSi/L
        end
    end
    

    kb_ini = griddata(dat_lon(~isnan(dat_kb)),dat_lat(~isnan(dat_kb)),dat_kb(~isnan(dat_kb)),lon,lat,'natural');
    kb_ini(mask==0) = 0;

    syn_ini = griddata(dat_lon(~isnan(dat_syn)),dat_lat(~isnan(dat_syn)),dat_syn(~isnan(dat_syn)),lon,lat,'natural');
    syn_ini(mask==0) = 0;

    sal_ini = griddata(dat_lon(~isnan(dat_sal)),dat_lat(~isnan(dat_sal)),dat_sal(~isnan(dat_sal)),lon,lat,'natural');
    sal_ini(mask==0) = 0;

    temp_ini = griddata(dat_lon(~isnan(dat_temp)),dat_lat(~isnan(dat_temp)),dat_temp(~isnan(dat_temp)),lon,lat,'natural');
    temp_ini(mask==0) = 0;

    dpo4_ini = griddata(dat_lon(~isnan(dat_dpo4)),dat_lat(~isnan(dat_dpo4)),dat_dpo4(~isnan(dat_dpo4)),lon,lat,'natural');
    dpo4_ini(mask==0) = 0;

    dnh4_ini = griddata(dat_lon(~isnan(dat_dnh4)),dat_lat(~isnan(dat_dnh4)),dat_dnh4(~isnan(dat_dnh4)),lon,lat,'natural');
    dnh4_ini(mask==0) = 0;

    dnox_ini = griddata(dat_lon(~isnan(dat_dnox)),dat_lat(~isnan(dat_dnox)),dat_dnox(~isnan(dat_dnox)),lon,lat,'natural');
    dnox_ini(mask==0) = 0;

    don_ini = griddata(dat_lon(~isnan(dat_don)),dat_lat(~isnan(dat_don)),dat_don(~isnan(dat_don)),lon,lat,'natural');
    don_ini(mask==0) = 0;

    dop_ini = griddata(dat_lon(~isnan(dat_dop)),dat_lat(~isnan(dat_dop)),dat_dop(~isnan(dat_dop)),lon,lat,'natural');
    dop_ini(mask==0) = 0;

    tn_ini = griddata(dat_lon(~isnan(dat_tn)),dat_lat(~isnan(dat_tn)),dat_tn(~isnan(dat_tn)),lon,lat,'natural');
    tn_ini(mask==0) = 0;

    tp_ini = griddata(dat_lon(~isnan(dat_tp)),dat_lat(~isnan(dat_tp)),dat_tp(~isnan(dat_tp)),lon,lat,'natural');
    tp_ini(mask==0) = 0;

    sio2_ini = griddata(dat_lon(~isnan(dat_sio2)),dat_lat(~isnan(dat_sio2)),dat_sio2(~isnan(dat_sio2)),lon,lat,'natural');
    sio2_ini(mask==0) = 0;

    figure;
    m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
    m_contourf(lon,lat,log10(syn_ini),100,'linestyle','none');
    hold on;
    m_gshhs_h('patch',[.6 .6 .6]);
    hold on;
    m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
    colorbar;
    set(gcf,'color', [1 1 1]);
    xlabel('Longitude','fontsize',9);
    ylabel('Latitude','fontsize',9);
end