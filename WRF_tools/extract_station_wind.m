function [u10,v10,slp,time] = extract_station_wind(res_path,slon,slat,domain_id)
%   Detailed explanation goes here

addpath(path,'/home/ychen/West_florida/Self_functions');


kk=0;

cmap = colormap(turbo);
nrep = 30;
nskip = 30;
cmap2(:,1) = interp1([1,nrep],[1,cmap(nrep,1)],1:nrep);
cmap2(:,2) = interp1([1,nrep],[1,cmap(nrep,2)],1:nrep);
cmap2(:,3) = interp1([1,nrep],[1,cmap(nrep,3)],1:nrep);
cmap = [cmap2;cmap((nrep+1):end,:)];
cmap = cmap(1:end-nskip,:);


for i = 1:1000

    date = datenum(2022,9,27)+i/24;
    fn = [res_path,'/wrfout_d',sprintf('%02d',domain_id),'_',datestr(date,'yyyy-mm-dd_HH:MM:SS')];

    if(exist(fn,'file'))
        if(i==1)
            lon = ncread(fn,'XLONG');
            lat = ncread(fn,'XLAT');
            [s_i,s_j] = find_ij(slon,slat,lon,lat,ones(size(lon)));
        end

        var1 = ncread(fn,'U10');
        var2 = ncread(fn,'V10');
        var3 = ncread(fn,'PSFC');

        time(i) = date;
        u10(i) = double(var1(s_i,s_j));
        v10(i) = double(var2(s_i,s_j));
        slp(i) = double(var3(s_i,s_j));
    else
        break;
    end
end
end