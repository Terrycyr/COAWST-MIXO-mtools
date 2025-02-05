clear all;

grd_name = 'cpb_grd_80x120_wf32.nc';

lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');
gn = struct('lon_rho',lon);

out_ij = [];
fid = fopen("ps_input",'rt+');
flag = 0;
k = 0;

while(~feof(fid))
    dat = fgetl(fid);
    if(contains(dat,'Location'))
        flag = 1;
        continue;
    end

    if(contains(dat,'-99'))
        break;
    end

    if(flag==1)
        k = k+1;
        tmp = str2num(dat(1:16));
        st_name{k} = strtrim(dat(17:end));
        out_ij(k,1) = tmp(2);
        out_ij(k,2) = tmp(3);
        fgetl(fid);
    end
end


figure;
contourf(lon,lat,mask);
hold on;
for i=1:length(st_name)
    scatter(lon(out_ij(i,1),out_ij(i,2)),lat(out_ij(i,1),out_ij(i,2)),40,'green','filled');
    hold on;
    text(lon(out_ij(i,1),out_ij(i,2)),lat(out_ij(i,1),out_ij(i,2)),st_name{i},'color','r','fontweight','bold');
end