% By Yuren Chen, 2021-04-14, Version 1.0
% Program to download ghrsst data and save to nc file by day
clear 

%Download time range, format: yyyymmdd
for year = 2001:2020
%year = 2001;
startdate = year*10000+0101
enddate = year*10000+1231

down_path = ['./',num2str(year),'/'];
if(~exist(down_path))
    mkdir(down_path);
end

%fn=['https://thredds.jpl.nasa.gov/thredds/dodsC/OceanTemperature/MUR-JPL-L4-GLOB-v4.1.nc'];
fn = ['https://thredds.jpl.nasa.gov/thredds/dodsC/OceanTemperature/MW_OI-REMSS-L4-GLOB-v5.0.nc']
try
    lon=double(ncread(fn,'lon'));
    lat=double(ncread(fn,'lat'));
    fprintf('lon & lat downloaded')
    fprintf('\n')
catch
    fprintf('Connection Failed')
    fprintf('\n')
end

tmp = find(lon>-90);
start(1) = tmp(1);
tmp = find(lat>20);
start(2) = tmp(1);

tmp = find(lon<-78);
count(1) = tmp(end) - start(1) +1;
tmp = find(lat<32);
count(2) = tmp(end) - start(2) +1;

time=double(ncread(fn,'time'));
nn=size(time,1)
tt=datenum(1981,1,1);

for i=1:nn;
    ts(i,:)=datestr(time(i)/3600/24+tt,'yyyymmdd');
    tnum(i) = time(i)/3600/24+tt;
    tn(i) = str2num(datestr(time(i)/24+tt,'yyyymmdd'));
end

tmp = 0;
k=0;
while(tmp(end)<nn)
    k=k+1;
    tp(k) = tmp(end)+1;
    tmp = find(tn == tn(tmp(end)+1));
end 
% %--------------------------------------------------------------------------

[xii,yii]=meshgrid(lon(start(1):start(1)+count(1)-1,1),...
                   lat(start(2):start(2)+count(2)-1,1));
%xii=xii';
%yii=yii';

xii = lon(start(1):start(1)+count(1)-1,1);
yii = lat(start(2):start(2)+count(2)-1,1);

xii=single(xii);
yii=single(yii);

for i=1:length(tp)
    if(str2num(ts(tp(i),:))>=startdate&&str2num(ts(tp(i),:))<=enddate)
        start(3) = tp(i);
        year_s=ts(start(3),1:4);
        ts(start(3),:)      
        if(i<length(tp))
            count(3) = tp(i+1)-tp(i);
        else
            count(3) = length(tn)-tp(i)+1;
        end
        flag = 0;
        while(flag==0)
        try
            sst=ncread(fn,'analysed_sst',start,count)-273.15;
            fprintf('sst downloaded')
            fprintf('\n')

            %asst=ncread(fn,'sst_anomaly',start,count);
            %fprintf('asst downloaded')
            %fprintf('\n')

            esst=ncread(fn,'analysis_error',start,count);
            fprintf('esst downloaded')
            fprintf('\n')

            mask=ncread(fn,'mask',start,count);
            fprintf('mask downloaded')
            fprintf('\n')
            
            flag = 1;
        catch
            flag = 0;
            "Error! Retrying..."
        end
        end

        %------------------------------------------------------------------------
        fn_nc=['./',down_path,'/GHRSST_MW_OI_REMSS_',ts(start(3),:),'.nc'];

        nccreate(fn_nc,'lon','Dimensions',{'lon', count(1)},'datatype','single');
        nccreate(fn_nc,'lat','Dimensions',{'lat', count(2)},'datatype','single');
        nccreate(fn_nc,'time','Dimensions',{'time', count(3) },'datatype','single');
        nccreate(fn_nc,'analysed_sst','Dimensions',{'lon', count(1), 'lat', count(2), 'time', count(3)},'datatype','single');
        %nccreate(fn_nc,'sst_anomaly','Dimensions',{'lon', count(1), 'lat', count(2), 'time', count(3)},'datatype','single');
        nccreate(fn_nc,'analysis_error','Dimensions',{'lon', count(1), 'lat', count(2), 'time', count(3)},'datatype','single');
        nccreate(fn_nc,'mask','Dimensions',{'lon', count(1), 'lat', count(2), 'time', count(3)},'datatype','single');

        sst=single(sst);
        %asst=single(asst);
        esst=single(esst);
        mask=single(mask);

        t = single (tnum(start(3):(start(3)+count(3)-1)));
        ncwrite(fn_nc,'lon',xii);
        ncwrite(fn_nc,'lat',yii);
        ncwrite(fn_nc,'time',t);
        ncwrite(fn_nc,'analysed_sst',sst);
        %ncwrite(fn_nc,'sst_anomaly',asst);
        ncwrite(fn_nc,'analysis_error',esst);
        ncwrite(fn_nc,'mask',mask);
    end
end
clear all
end
