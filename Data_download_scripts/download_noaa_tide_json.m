clearvars;

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
    format = 'json';

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

    json = webread(url);

    try 
        station_name= json.metadata.name;
        station_id= json.metadata.id;
        station_lat= json.metadata.lat;
        station_lon= json.metadata.lon;
        k=k+1;
        disp(station_name)
    catch
        disp(['No data for selected time at station ',station(i,:),'!'])
        continue;  
    end

    station_name= json.metadata.name;
    station_id= json.metadata.id;
    station_lat= json.metadata.lat;
    station_lon= json.metadata.lon;

    dat_raw2 = struct2table(json.data);

    id(k) = str2num(station_id);
    lat(k) = str2num(station_lat);
    lon(k) = str2num(station_lon);
    name{k} = station_name;
    time{k} = datenum(cell2mat(table2array(dat_raw2(:,"t"))));
    time_vec{k} = datevec(time{k});
    water_level{k} = table2array(dat_raw2(:,"v"));
end

save(save_name,"id","lat","lon","name","time","time_vec","water_level");