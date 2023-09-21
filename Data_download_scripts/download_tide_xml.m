clear all;

product = 'hourly_height';
application = 'NOS.COOPS.TAC.WL';
begin_date = '20220101';
end_date = '20230101';
datum = 'MSL';
station = num2str(load("WF_stations.txt"));
time_zone = 'GMT';
units = 'metric';
save_name = 'NOAA_Tide_WFS_2022.mat';

k=0;
for i = 1:size(station,1)

    %Retrieve Metadata
    format = 'xml';

    url = ['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?' ...
        ,'product=',product ...
        ,'&application=',application ...
        ,'&begin_date=',begin_date ...
        ,'&end_date=',end_date ...
        ,'&datum=',datum ...
        ,'&station=',station(i,:) ...
        ,'&time_zone=',time_zone ...
        ,'&units=',units ...
        '&format=',format];

    xml = webread(url);
    fid = fopen('tmp.xml','wt+');
    fprintf(fid,'%s',xml);
    fclose(fid);
    
    dat_raw1 = parseXML('tmp.xml');
    [station_id,station_lat,station_lon,station_name] = dat_raw1.Children(2).Attributes.Value;

    delete('tmp.xml');

    %Retrieve dataset
    format = 'csv';

    url = ['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?' ...
        ,'product=',product ...
        ,'&application=',application ...
        ,'&begin_date=',begin_date ...
        ,'&end_date=',end_date ...
        ,'&datum=',datum ...
        ,'&station=',station(i,:) ...
        ,'&time_zone=',time_zone ...
        ,'&units=',units ...
        '&format=',format];

    dat_raw2 = webread(url);

    if(size(dat_raw2,1)==0)
        disp(['No data for selected time at station ',station_name,' ',station_id,'!'])
        continue;
    else
        k=k+1;
        disp(station_name)
    end
    id(k) = str2num(station_id);
    slat(k) = str2num(station_lat);
    slon(k) = str2num(station_lon);
    sname{k} = station_name;
    obs_time{k} = datenum(table2array(dat_raw2(:,"DateTime")));
    obs_time_vec{k} = datevec(obs_time{k});
    obs_water_level{k} = table2array(dat_raw2(:,"WaterLevel"));
end

save(save_name,"id","slat","slon","sname","obs_time","obs_time_vec","obs_water_level");