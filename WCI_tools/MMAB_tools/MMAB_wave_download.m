% Programmed by terry 2017/10/6.
% To download files from the ftp of noaa.
% Without any change, the MMAB data during the whole 2010 will be downloaded.

clear all
year = 2017;
ncep = ftp('polar.ncep.noaa.gov');
for i=1:12
    i
    if(i<10)
        mydate = [num2str(year),'0',num2str(i)];
    else
        mydate = [num2str(year),num2str(i)];
    end
    cd(ncep,['/pub/history/waves/multi_1/', mydate, '/gribs']);
    mget(ncep,['multi_1.glo_30m.dp.',mydate,'.grb2']);
    mget(ncep,['multi_1.glo_30m.hs.',mydate,'.grb2']);
    mget(ncep,['multi_1.glo_30m.tp.',mydate,'.grb2']);
    %mget(ncep,['multi_1.glo_30m.wind.',mydate,'.grb2']);
end