% By Weicong Cheng, 2016-09-22, Version 1.0
% Program to download HYCOM data and save to nc file by day

% Modified by Weicong Cheng, 2016-10-29, Version 2.0
% Download count(3) layers' data in one time

% Modified by Weicong Cheng, 2017-09-03, Version 2.1
% Read the time from the data

% Modified by Weicong Cheng, 2017-09-04, Version 3.0
% For the data from 2013/01/01 to present

% Modified by Yuren Chen, 2020-08-27, Version 4.0
% Modified by Yuren Chen, 2022-12-19, Version 5.0, fixed bugs of missing time records in the first file.
% Modified by Yuren Chen, 2023-03-24, Version 5.1, added basic quality checking codes.

clear

myName= mfilename;
%Download time range, format: yyyymmdd
startdate = 20220127
enddate = 20230101
lon_range = [-90 -78];
lat_range = [20 32];
fday_treat_flag = 0;

%url = 'http:///tds.hycom.org/thredds/dodsC/GLBv0.08/expt_53.X/data/';
%url = 'https://tds.hycom.org/thredds/dodsC/GLBv0.08/expt_93.0';
url = 'https://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0';

%--------------------- NO CHANGE BELOW-------------------------------------
year = floor(startdate/10000);
down_path = ['./',num2str(year),'/'];
if(~exist(down_path))
    mkdir(down_path);
end

%fn=[url, num2str(year)];
fn = url;

try
    lon=ncread(fn,'lon');
    lat=ncread(fn,'lat');
    fprintf('lon & lat downloaded')
    fprintf('\n')
catch
    fprintf('Connection Failed!')
    fprintf('\n')
    return
end

if(~sum(lon<0))
    if(lon_range(1)<0)
        lon_range(1) = 360+lon_range(1);
    end
    if(lon_range(2)<0)
        lon_range(2) = 360+lon_range(2);
    end
    lon_range = sort(lon_range);
end

tmp = find(lon>lon_range(1));
start(1) = tmp(1);
tmp = find(lat>lat_range(1));
start(2) = tmp(1);

tmp = find(lon<lon_range(2));
count(1) = tmp(end) - start(1) +1;
tmp = find(lat<lat_range(2));
count(2) = tmp(end) - start(2) +1;

%Total Layers
dep=ncread(fn,'depth');
fprintf('depth downloaded')
fprintf('\n')

start(3)=1;
count(3)=length(dep);

time=ncread(fn,'time');
nn=size(time,1);
tt=datenum(2000,1,1);

for i=1:nn
    ts(i,:)=datestr(time(i)/24+tt,'yyyymmdd');
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

% [xii,yii]=meshgrid(lon(start(1):start(1)+count(1)-1,1),...
%                    lat(start(2):start(2)+count(2)-1,1));
% xii=xii';
% yii=yii';
xii=lon(start(1):start(1)+count(1)-1,1);
yii=lat(start(2):start(2)+count(2)-1,1);

xii=single(xii);
yii=single(yii);
dep=single(dep);

for i=1:length(tp)
    if(str2num(ts(tp(i),:))>=startdate&&str2num(ts(tp(i),:))<=enddate)
        tmp = datevec(datenum(2000,1,1)+time(tp(i))/24);
        dat_hour = tmp(4);
        start(4) = tp(i);
        year_s=ts(start(4),1:4);
        ts(start(4),:)
        if(i<length(tp))
            count(4) = tp(i+1)-tp(i);
        else
            count(4) = length(tn)-tp(i)+1;
        end

        t = single (time(start(4):(start(4)+count(4)-1)));
        % special treat for the first day !!!
        if(i==1&&dat_hour>0&&fday_treat_flag==1)
            t_b = datevec(datenum(num2str(tn(i)),'yyyymmdd')-1);
            year_b = t_b(1);
            if(year_b~=year)
                fn_b = [url, num2str(year_b)];
                time_b = ncread(fn_b,'time');
            else
                fn_b = fn;
                time_b = time;
            end
            nn=size(time_b,1);
            for ii=1:nn
                ts_b(ii,:)=datestr(time_b(ii)/24+tt,'yyyymmdd');
                tn_b(ii) = str2num(datestr(time_b(ii)/24+tt,'yyyymmdd'));
            end
            tmp = find(tn_b == tn(i));
            tp_b = tmp(1);
            count_b = length(tmp);

            t_b = time_b(tp_b:(tp_b+count_b-1));

            tem_b=ncread(fn_b,'water_temp',...
                [start(1) start(2) start(3) tp_b],[count(1) count(2) count(3) count_b]);
            fprintf('temperature downloaded')
            fprintf('\n')

            u_b=ncread(fn_b,'water_u', ...
                [start(1) start(2) start(3) tp_b],[count(1) count(2) count(3) count_b]);
            fprintf('u downloaded')
            fprintf('\n')

            v_b=ncread(fn_b,'water_v', ...
                [start(1) start(2) start(3) tp_b],[count(1) count(2) count(3) count_b]);
            fprintf('v downloaded')
            fprintf('\n')

            sal_b=ncread(fn_b,'salinity', ...
                [start(1) start(2) start(3) tp_b],[count(1) count(2) count(3) count_b]);
            fprintf('salinity downloaded')
            fprintf('\n')


            el_b=ncread(fn_b,'surf_el',[start(1) start(2) tp_b],[count(1) count(2) count_b]);
            fprintf('surf_el downloaded')
            fprintf('\n')
        end
        flag = 0;
        count2 = 0;
        while(flag==0)
             try
                tem=ncread(fn,'water_temp',start,count);
                fprintf('temperature downloaded')
                fprintf('\n')
                flag = 1;

                tmp = tem(~isnan(tem));
                maxval = max(abs(tmp));
                if(maxval>100)
                    count2 = count2+1;
                    if(count2>5)
                        'Bad Records detected! Skipped!'
                        flag = 1;
                    else
                        flag = 0;
                    end
                end
            catch
                'Errors in downloading temperature!'
                flag = 0;
            end
        end

        if(count2>5)
            break;
        end

        flag = 0;
        count2 = 0;
        while(flag==0)
             try
                u=ncread(fn,'water_u',start,count);
                fprintf('u downloaded')
                fprintf('\n')
                flag = 1;

                tmp = u(~isnan(u));
                maxval = max(abs(tmp));
                if(maxval>100)
                    count2 = count2+1;
                    if(count2>5)
                        'Bad Records detected! Skipped!'
                        flag = 1;
                    else
                        flag = 0;
                    end
                end
            catch
                'Errors in downloading u!'
                flag = 0;
            end
        end

        if(count2>5)
            break;
        end

        flag = 0;
        count2 = 0;
        while(flag==0)
            try
                v=ncread(fn,'water_v',start,count);
                fprintf('v downloaded')
                fprintf('\n')
                flag = 1;

                tmp = v(~isnan(v));
                maxval = max(abs(tmp));
                if(maxval>100)
                    count2 = count2+1;
                    if(count2>5)
                        'Bad Records detected! Skipped!'
                        flag = 1;
                    else
                        flag = 0;
                    end
                end
            catch
                'Errors in downloading v!'
                flag = 0;
            end
        end

        if(count2>5)
            break;
        end

        flag = 0;
        count2 = 0;
        while(flag==0)
            try
                sal=ncread(fn,'salinity',start,count);
                fprintf('salinity downloaded')
                fprintf('\n')
                flag = 1;

                tmp = sal(~isnan(sal));
                maxval = max(abs(tmp));
                if(maxval>100)
                    count2 = count2+1;
                    if(count2>5)
                        'Bad Records detected! Skipped!'
                        flag = 1;
                    else
                        flag = 0;
                    end
                end
            catch
                'Errors in downloading salinity!'
                flag = 0;
            end
        end

        if(count2>5)
            break;
        end

        flag = 0;
        count2 = 0;
        while(flag==0)
            try
                el=ncread(fn,'surf_el',[start(1) start(2) start(4)],[count(1) count(2) count(4)]);
                fprintf('surf_el downloaded')
                fprintf('\n')
                flag = 1;

                tmp = el(~isnan(el));
                maxval = max(abs(tmp));
                if(maxval>100)
                    count2 = count2+1;
                    if(count2>5)
                        'Bad Records detected! Skipped!'
                        flag = 1;
                    else
                        flag = 0;
                    end
                end
            catch
                'Errors in downloading el!'
                flag = 0;
            end
        end 

        if(count2>5)
            break;
        end

        if(exist('t_b','var'))
            d1 = length(t_b);
            d2 = length(t);

            t2(1:d1) = t_b;
            t2([1:d1]+d1) = t;

            tem2(:,:,:,1:d1) = tem_b;
            tem2(:,:,:,[1:d2]+d1) = tem;

            u2(:,:,:,1:d1) = u_b;
            u2(:,:,:,[1:d2]+d1) = u;

            v2(:,:,:,1:d1) = v_b;
            v2(:,:,:,[1:d2]+d1) = v;

            sal2(:,:,:,1:d1) = sal_b;
            sal2(:,:,:,[1:d2]+d1) = sal;

            el2(:,:,1:d1) = el_b;
            el2(:,:,[1:d2]+d1) = el;

            t = t2;
            tem = tem2;
            u = u2;
            v=v2;
            sal = sal2;
            el = el2;

            count(4) = count(4)+count_b;
            clear t_b tem_b u_b v_b sal_b el_b
        end
        %------------------------------------------------------------------------
        fn_nc=['./',down_path,'/HYCOM_',ts(start(4),:),'.nc'];

        ncid=netcdf.create(fn_nc,'clobber');
        disp(' ## Defining Global Attributes...')
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by ',myName,' on ' datestr(now)]);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'type', ...
            ['HYCOM file from',url]);
        netcdf.close(ncid);

        nccreate(fn_nc,'lon','Dimensions',{'lon',count(1)},'datatype','single');
        nccreate(fn_nc,'lat','Dimensions',{'lat',count(2)},'datatype','single');
        nccreate(fn_nc,'time','Dimensions',{'time',count(4)},'datatype','single');
        nccreate(fn_nc,'depth','Dimensions',{'layer',count(3)},'datatype','single');
        nccreate(fn_nc,'temperature','Dimensions',{'lon',count(1),'lat',count(2), 'layer',count(3),'time',count(4)},'datatype','single');
        nccreate(fn_nc,'u','Dimensions',{'lon',count(1),'lat',count(2), 'layer',count(3),'time',count(4)},'datatype','single');
        nccreate(fn_nc,'v','Dimensions',{'lon',count(1),'lat',count(2), 'layer',count(3),'time',count(4)},'datatype','single');
        nccreate(fn_nc,'salinity','Dimensions',{'lon',count(1),'lat',count(2), 'layer',count(3),'time',count(4)},'datatype','single');
        nccreate(fn_nc,'surf_el','Dimensions',{'lon',count(1),'lat',count(2),'time',count(4)},'datatype','single');


        tem=single(tem);
        u=single(u);
        v=single(v);
        sal=single(sal);
        el=single(el);

        xii(xii>180) = xii(xii>180)-360;

        ncwrite(fn_nc,'lon',xii);
        ncwrite(fn_nc,'lat',yii);
        ncwrite(fn_nc,'time',t);
        ncwrite(fn_nc,'depth',dep);
        ncwrite(fn_nc,'temperature',tem);
        ncwrite(fn_nc,'u',u);
        ncwrite(fn_nc,'v',v);
        ncwrite(fn_nc,'salinity',sal);
        ncwrite(fn_nc,'surf_el',el);
    end
end
