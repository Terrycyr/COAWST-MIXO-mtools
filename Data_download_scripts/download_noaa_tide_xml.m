clear all;

product = 'hourly_height';
application = 'NOS.COOPS.TAC.WL';
begin_date = '20030101';
end_date = '20040101';
datum = 'MSL';
station = num2str(load("noaa_station_32.txt"));
time_zone = 'GMT';
units = 'metric';
save_name = 'NOAA_Tide_CPB.mat';
k=0;
for i = 1:size(station,1)

    %Retrieve Metadata
    format = 'xml';

    url = ['https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?' ...
        ,'product=',product ...
        ,'&application=',application ...
        ,'&begin_date=',begin_date ...
        ,'&end_date=',begin_date ...
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
    lat(k) = str2num(station_lat);
    lon(k) = str2num(station_lon);
    name{k} = station_name;
    time{k} = datenum(table2array(dat_raw2(:,"DateTime")));
    time_vec{k} = datevec(time{k});
    water_level{k} = table2array(dat_raw2(:,"WaterLevel"));
end

save(save_name,"id","lat","lon","name","time","time_vec","water_level");