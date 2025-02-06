clear all;close all;
addpath(path,'/home/ychen/West_florida/Self_functions');

station = {'C10','C12','C13'};

u_path='/home/ychen/West_florida/2022_IAN_WCI/WFS_2022_Uwind1.nc';
v_path='/home/ychen/West_florida/2022_IAN_WCI/WFS_2022_Vwind1.nc';
slp_path='/home/ychen/West_florida/2022_IAN_WCI/WFS_2022_Pair1.nc';
grd_path='/home/ychen/West_florida/2022_IAN_WCI/ROMS_WFS_new.nc';

figure('Units','pixels','Position',[100 100 1000 550]);
tiledlayout(3,2,'TileSpacing','compact');

for i_num=1:3
  

    %extract data from buoy
    [~,~,dat1] = xlsread(['./Ian_bouy/buoy_',station{i_num},'_wind_spd.xlsx']);
    [~,~,dat2] = xlsread(['./Ian_bouy/buoy_',station{i_num},'_wind_dir.xlsx']);

    loc_raw = dat1{2,1};
    loc1 = strfind(loc_raw,'(');
    loc2 = strfind(loc_raw,'N');
    loc3 = strfind(loc_raw,'E');
    slat = str2num(loc_raw(loc1+1:loc2-1));
    slon = str2num(loc_raw(loc2+1:loc3-1));

    k = 0;

    for i=1:size(dat1,1)
        if(isnumeric(dat1{i,1}))
            k=k+1;
            if(strcmp(dat1{i,3},'good'))
                buoy_time(k) = xlstime2date(dat1{i,1});
                buoy_spd0(k) = dat1{i,2};
            else
                buoy_time(k) = xlstime2date(dat1{i,1});
                buoy_spd0(k) = NaN;
            end
        end
    end

    k2 = 0;

    for i=1:size(dat2,1)
        if(isnumeric(dat2{i,1}))
            k2=k2+1;
            if(strcmp(dat2{i,3},'good'))
                buoy_time2(k2) = xlstime2date(dat2{i,1});
                buoy_dir0(k2) = dat2{i,2};
            else
                buoy_time2(k2) = xlstime2date(dat2{i,1});
                buoy_dir0(k2) = NaN;
            end
        end
    end
    %end
    

    %[u10,v10,slp,time] = extract_station_wind(result_path,slon,slat,domain_id);
    [u10,v10,slp,time] = extract_station_wind_2(u_path,v_path,slp_path,grd_path,slon,slat);
    time = time+datenum(2022,1,1);

    spd = sqrt(u10.^2+v10.^2);
    dir = 90-(atan2d(-1*v10,-1*u10));
    dir(dir<0) = dir(dir<0)+360;

%     [u10_GFS,v10_GFS,slp_GFS,time_GFS] = extract_gfs_wind(result_path_GFS,slon,slat);
%     spd_GFS0 = sqrt(u10_GFS.^2+v10_GFS.^2);

    % Hsu, S. A., Eric A. Meindl, and David B. Gilhousen, 1994: Determining the Power-Law Wind-Profile Exponent under Near-Neutral Stability Conditions at Sea, Applied Meteorology, Vol. 33, No. 6, June 1994.
    buoy_spd_10 = buoy_spd0*(10/3.1)^0.11;
    buoy_spd = interp1(buoy_time,buoy_spd_10,time);

    buoy_dir = interp1(buoy_time2,buoy_dir0,time);

    x1 = spd;
    x2 = buoy_spd;
    x1(isnan(x2)) = [];
    x2(isnan(x2)) = [];
    tmp = corrcoef(x1,x2);
    CC(i_num) = tmp(1,2);
    RMSE(i_num) = cal_rmse(x1,x2);

    nexttile;
    plot(time,spd,'color',[0.25 0.25 0.25],'LineWidth',2);
    hold on;
    scatter(buoy_time,buoy_spd_10,30,[0.8500 0.3250 0.0980],'filled');
    hold on;
    alpha(0.4);
    text('Units','normalized','Position',[0.03 0.9],'String',station{i_num},'FontSize',12);
    hold on;
    text('Units','normalized','Position',[0.03 0.68],'String',strcat('RMSE=',num2str(RMSE(i_num),2)),'FontSize',10);
    hold on;
    text('Units','normalized','Position',[0.03 0.55],'String',strcat('CC=',num2str(CC(i_num),2)),'FontSize',10);
    xlim([datenum(2022,9,27) datenum(2022,10,2)]);
    if(i_num==3)
        set(gca,'xtick',[datenum(2022,9,27):1:datenum(2022,10,2)],...
            'xticklabel',datestr([datenum(2022,9,27):1:datenum(2022,10,2)],'mm/dd'));
        xlabel('Time (GMT)');
    else
        set(gca,'xtick',[datenum(2022,9,27):12/24:datenum(2022,10,2)],'XTickLabel',[]);
    end
    
    set(gca,'fontsize',11);

    axes_pos = get(gca,'position');

    annotation('textbox',[axes_pos(1)-0.07,axes_pos(2)+axes_pos(4)+0.02,0,0],...
        'String',['(',char('a'+(i_num-1)*2),')'],'fontsize',12);

    ylim([0 50]);
    ylabel('Speed m/s');

    if(i_num==1)
        hl=legend('WRF-GFS','Bouy');
        hl.NumColumns = 1;
        hl.Box = "off";
        hl.FontSize = 12;
    end
    set(gcf,'color','w');

    nexttile;
    plot(time,dir,'color',[0.25 0.25 0.25],'LineWidth',2);
    hold on;
    scatter(buoy_time2,buoy_dir0,30,[0.8500 0.3250 0.0980],'filled');
    hold on;
    alpha(0.4);
    text('Units','normalized','Position',[0.03 0.9],'String',station{i_num},'FontSize',12);
    hold on;
    xlim([datenum(2022,9,27) datenum(2022,10,2)]);
    if(i_num==3)
        set(gca,'xtick',[datenum(2022,9,27):1:datenum(2022,10,2)],...
            'xticklabel',datestr([datenum(2022,9,27):1:datenum(2022,10,2)],'mm/dd'));
        xlabel('Time (GMT)');
    else
        set(gca,'xtick',[datenum(2022,9,27):12/24:datenum(2022,10,2)],'XTickLabel',[]);
    end

    set(gca,'fontsize',11);

    axes_pos = get(gca,'position');

    annotation('textbox',[axes_pos(1)-0.07,axes_pos(2)+axes_pos(4)+0.02,0,0],...
        'String',['(',char('b'+(i_num-1)*2),')'],'fontsize',12);

    ylim([0 360]);
    ylabel('Direction ^{\circ}');

    set(gcf,'color','w');
    set(gca,'ytick',0:120:360);

    clear buoy_time buoy_time2 buoy_spd0 buoy_dir0
end

exportgraphics(gcf,'figxxx_Ian_verify_merged_spd_dir.png','Resolution',600);
