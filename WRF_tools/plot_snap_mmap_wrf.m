clear all; close all;

addpath(path,'/home/ychen/matlab_tools/m_map');

filename = 'Ian_D2.gif';

load('/home/ychen/WRF-GFS-Hurricane/HURDAT2/IAN2022.mat');

ian_t = IAN2022{2};
ian_lon = IAN2022{5};
ian_lat = IAN2022{6};

kk=0;

cmap = colormap(turbo);
nrep = 30;
nskip = 30;
cmap2(:,1) = interp1([1,nrep],[1,cmap(nrep,1)],1:nrep);
cmap2(:,2) = interp1([1,nrep],[1,cmap(nrep,2)],1:nrep);
cmap2(:,3) = interp1([1,nrep],[1,cmap(nrep,3)],1:nrep);
cmap = [cmap2;cmap((nrep+1):end,:)];
cmap = cmap(1:end-nskip,:);


for i = 1:576

    date = datenum(2022,9,27)+i/24;
    fn = ['wrfout_d02_',datestr(date,'yyyy-mm-dd_HH:MM:SS')];

    if(i==1)
        lon = ncread(fn,'XLONG');
        lat = ncread(fn,'XLAT');
    end

    var1 = ncread(fn,'U10');
    var2 = ncread(fn,'V10');
    var = sqrt(var1.^2+var2.^2);

    kk = kk+1;
    
    figure(2);

    colormap(cmap);

    m_proj('Mercator','lat',[24 30],'long',[-87 -80]);

    hp = m_contourf(lon,lat,var,100,'linestyle','none');
    colorbar;
    m_gshhs_h('patch',[.8 .8 .8],'linestyle','none');
    m_grid('linestyle','none','fontsize',11);
    title(datestr(date,'yyyy-mmm-dd HH:MM'));
    clim([0 40]);
    xlabel('Longitude');
    ylabel('Latitude');
    set(gcf,'color','w');
    hold on;
    m_plot(ian_lon,ian_lat,'color',[.5 .5 .5],'linewidth',2);
    hold on;
    m_scatter(ian_lon,ian_lat,30,'k','filled');
    hold on;
    ian_now_x = interp1(datenum(ian_t),ian_lon,date);
    ian_now_y = interp1(datenum(ian_t),ian_lat,date);
    m_scatter(ian_now_x,ian_now_y,80,'b','filled');
    hold on;
    m_text(ian_lon(1:2:end)-2,ian_lat(1:2:end),datestr(ian_t(1:2:end),'mmm-dd HH:MM'));
    hold off;

    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if kk == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end