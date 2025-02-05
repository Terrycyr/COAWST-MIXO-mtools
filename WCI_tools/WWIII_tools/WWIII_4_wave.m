% Programmed by terry 2017/10/6.
% To combine and process files from the ftp of noaa into Tar
clear all
year = 2022;
b_points = load('WW3_boudary_point.dat');
[r_b,c_b] = size(b_points);

%list for included WW3 files
hycom_dir = ['./',num2str(year),'/'];
fname = dir([hycom_dir 'WW3*.nc']);
for i = 1:size(fname,1)
    filename{i,1} = strcat(hycom_dir,fname(i).name);
end

for n_b=1:r_b 
    n_b
    station_y = b_points(n_b,2);
    station_x = b_points(n_b,1);
    station_dd = 2;
    fid = fopen(['./WW3_WAVE_',num2str(year),'_',num2str(n_b),'.bnd'],'wt+');
    fprintf(fid,'TPAR\n');
    for ifile = 1:size(filename,1)
        fn = filename{ifile};
        time = double(ncread(fn,'time'));
        [rt,ct] = size(time);
        for t = 1:rt
            dp = double(ncread(fn,'Tdir',[1 1 1 t],[Inf Inf Inf 1]));
            tp = double(ncread(fn,'Tper',[1 1 1 t],[Inf Inf Inf 1]));
            hs = double(ncread(fn,'Thgt',[1 1 1 t],[Inf Inf Inf 1]));
            station_dp = dp(station_y,station_x);
            station_tp = tp(station_y,station_x);
            station_hs = hs(station_y,station_x);

            station_time = datestr(time(t)/24+datenum(2013,1,1),'yyyymmdd.HHMMSS');

            if(~(isnan(station_hs)||isnan(station_tp)||isnan(station_dp)))
                fprintf(fid,'%s%5.1f%5.1f%7.1f%5.1f\n',station_time,station_hs...
                    ,station_tp,station_dp,station_dd);
            end
        end
    end
    fclose(fid);      
end