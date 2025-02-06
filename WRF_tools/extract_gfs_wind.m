function [u10,v10,slp,time] = extract_gfs_wind(path,slon,slat)
%   Detailed explanation goes here

addpath(path,'/home/ychen/West_florida/Self_functions');

for i = 1:1000

    date = datenum(2022,9,27)+(i-1)*6/24;
    fn = [path,'/gfs_4_',datestr(date,'yyyymmdd_HHMM'),'_000.nc'];

    if(exist(fn,'file'))
        if(i==1)
            lon0 = double(ncread(fn,'lon_0'));
            lon0(lon0>180) = lon0(lon0>180)-360;
            lat0 = double(ncread(fn,'lat_0'));
            [lat,lon] = meshgrid(lat0,lon0);
            [s_i,s_j] = find_ij(slon,slat,lon,lat,ones(size(lon)));
        end

        var1 = ncread(fn,'UGRD_P0_L103_GLL0',[1 1 1],[Inf Inf 1]);
        var2 = ncread(fn,'VGRD_P0_L103_GLL0',[1 1 1],[Inf Inf 1]);
        var3 = ncread(fn,'PRMSL_P0_L101_GLL0');

        time(i) = date;
        u10(i) = double(var1(s_i,s_j));
        v10(i) = double(var2(s_i,s_j));
        slp(i) = double(var3(s_i,s_j));
    else
        break;
    end
end
end