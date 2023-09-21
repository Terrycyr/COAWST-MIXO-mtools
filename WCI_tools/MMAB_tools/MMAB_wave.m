% Programmed by terry 2017/10/6.
% To combine and process files from the ftp of noaa into Tar
clear all
setup_nctoolbox
year = 2017;
b_points = load('MMAB_boudary_point.dat');
[r_b,c_b] = size(b_points);
for n_b=1:r_b 
    station_y = b_points(n_b,2);
    station_x = b_points(n_b,1);
    station_dd = 2;
    fid = fopen(['./MMAB_WAVE_',num2str(year),'_',num2str(n_b),'.bnd'],'wt+');
    fprintf(fid,'TPAR\n');
    for month = 7:9 %!!!!
        if(month<10)
            mydate = [num2str(year),'0',num2str(month)];
        else
            mydate = [num2str(year),num2str(month)];
        end
        nc_dp = ncgeodataset(['./',num2str(year),'/multi_1.glo_30m.dp.',mydate,'.grb2']);
        nc_tp = ncgeodataset(['./',num2str(year),'/multi_1.glo_30m.tp.',mydate,'.grb2']);
        nc_hs = ncgeodataset(['./',num2str(year),'/multi_1.glo_30m.hs.',mydate,'.grb2']);
        var_list_dp = nc_dp.variables;
        var_list_tp = nc_tp.variables;
        var_list_hs = nc_hs.variables;
        %lat = double(nc_dp{'lat'}(:));
        %lon = double(nc_dp{'lon'}(:));
        time = double(nc_dp{'time'}(:));
        [rt,ct] = size(time);
        for t = 1:rt-1
            dp = double(squeeze(nc_dp{var_list_dp{1}}(t,:,:)));
            tp = double(squeeze(nc_tp{var_list_tp{1}}(t,:,:)));
            hs = double(squeeze(nc_hs{var_list_hs{1}}(t,:,:)));
            station_dp = dp(station_y,station_x);
            station_tp = tp(station_y,station_x);
            station_hs = hs(station_y,station_x);
            if(mod(time(t),24)<10)
                station_time = [num2str(int32(year*10000+month*100+floor(1+(time(t)...
                    /24)))),'.0',num2str(mod(time(t),24)),'0000'];
            else
                station_time = [num2str(int32(year*10000+month*100+floor(1+(time(t)...
                    /24)))),'.',num2str(mod(time(t),24)),'0000'];
            end
            if(~(isnan(station_hs)||isnan(station_tp)||isnan(station_dp)))
                fprintf(fid,'%s%5.1f%5.1f%7.1f%5.1f\n',station_time,station_hs...
                    ,station_tp,station_dp,station_dd);
            end
        end
    end
    fclose(fid);      
end