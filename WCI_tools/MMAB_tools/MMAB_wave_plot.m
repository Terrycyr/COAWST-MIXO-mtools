% Programmed by terry 2017/10/6.
% To combine and process files from the ftp of noaa into Tar
clear all
setup_nctoolbox
year = 2010;
for month = 10:10
    if(month<10)
        mydate = [num2str(year),'0',num2str(month)]
    else
        mydate = [num2str(year),num2str(month)]
    end
    nc_dp = ncgeodataset(['./multi_1.glo_30m.dp.',mydate,'.grb2']);
    nc_tp = ncgeodataset(['./multi_1.glo_30m.tp.',mydate,'.grb2']);
    nc_hs = ncgeodataset(['./multi_1.glo_30m.hs.',mydate,'.grb2']);
    var_list_dp = nc_dp.variables;
    var_list_tp = nc_tp.variables;
    var_list_hs = nc_hs.variables;
    lat = double(nc_dp{'lat'}(:));
    lon = double(nc_dp{'lon'}(:));
    time = double(nc_dp{'time'}(:));
    [rt,ct] = size(time);
    for t = 1:1:rt-1
        dp = double(squeeze(nc_dp{var_list_dp{1}}(t,:,:)));
        tp = double(squeeze(nc_tp{var_list_tp{1}}(t,:,:)));
        hs = double(squeeze(nc_hs{var_list_hs{1}}(t,:,:)));
        figure(1);
        contourf(lon,lat,hs,0:0.5:5);
        colorbar;
        axis equal;
        axis([108 122 18 25]);
    end
end