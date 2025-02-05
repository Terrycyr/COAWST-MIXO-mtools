clear all

myName= mfilename;
%Download time range, format: yyyymmdd
startdate = 20230101
enddate = 20240101
lon_range = [-90 -78];
lat_range = [20 32];
fday_treat_flag = 0;

url = 'https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ww3_global/WaveWatch_III_Global_Wave_Model_best.ncd';

%--------------------- NO CHANGE BELOW-------------------------------------
year = floor(startdate/10000);
down_path = ['./',num2str(year),'/'];
if(~exist(down_path))
    mkdir(down_path);
end

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
dep=ncread(fn,'z');
fprintf('depth downloaded')
fprintf('\n')

start(3)=1;
count(3)=length(dep);

time=ncread(fn,'time');
nn=size(time,1);
tt=datenum(2013,1,1);

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
        tmp = datevec(datenum(2013,1,1)+time(tp(i))/24);
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

            hs_b=ncread(fn_b,'Thgt',...
                [start(1) start(2) start(3) tp_b],[count(1) count(2) count(3) count_b]);
            fprintf('Significant wave height downloaded')
            fprintf('\n')

            wper_b=ncread(fn_b,'Tper', ...
                [start(1) start(2) start(3) tp_b],[count(1) count(2) count(3) count_b]);
            fprintf('Wave period downloaded')
            fprintf('\n')

            wdir_b =ncread(fn_b,'Tdir', ...
                [start(1) start(2) start(3) tp_b],[count(1) count(2) count(3) count_b]);
            fprintf('Wave direction downloaded')
            fprintf('\n')
        end
        flag = 0;
        count2 = 0;
        while(flag==0)
             try
                hs=ncread(fn,'Thgt',start,count);
                fprintf('Significant wave height downloaded')
                fprintf('\n')
                flag = 1;

                tmp = hs(~isnan(hs));
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
                'Errors in downloading significant wave height!'
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
                wper=ncread(fn,'Tper',start,count);
                fprintf('Wave period downloaded')
                fprintf('\n')
                flag = 1;

                tmp = wper(~isnan(wper));
                maxval = max(abs(tmp));
                if(maxval>999)
                    count2 = count2+1;
                    if(count2>5)
                        'Bad Records detected! Skipped!'
                        flag = 1;
                    else
                        flag = 0;
                    end
                end
            catch
                'Errors in downloading wave period!'
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
                wdir=ncread(fn,'Tdir',start,count);
                fprintf('Wave direction downloaded')
                fprintf('\n')
                flag = 1;

                tmp = wdir(~isnan(wdir));
                maxval = max(abs(tmp));
                if(maxval>999)
                    count2 = count2+1;
                    if(count2>5)
                        'Bad Records detected! Skipped!'
                        flag = 1;
                    else
                        flag = 0;
                    end
                end
            catch
                'Errors in downloading wave direction!'
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

            hs2(:,:,:,1:d1) = hs_b;
            hs2(:,:,:,[1:d2]+d1) = hs;

            wper2(:,:,:,1:d1) = wper_b;
            wper2(:,:,:,[1:d2]+d1) = wper;

            wdir2(:,:,:,1:d1) = wdir_b;
            wdir2(:,:,:,[1:d2]+d1) = wdir;

            t = t2;
            hs = hs2;
            wper = wper2;

            count(4) = count(4)+count_b;
            clear t_b hs_b wper_b wdir_b
        end
        %------------------------------------------------------------------------
        fn_nc=['./',down_path,'/WW3_',ts(start(4),:),'.nc'];

        ncid=netcdf.create(fn_nc,'clobber');
        disp(' ## Defining Global Attributes...')
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by ',myName,' on ' datestr(now)]);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'type', ...
            ['WW3 file from',url]);
        netcdf.close(ncid);

        nccreate(fn_nc,'lon','Dimensions',{'lon',count(1)},'datatype','single');
        nccreate(fn_nc,'lat','Dimensions',{'lat',count(2)},'datatype','single');
        nccreate(fn_nc,'time','Dimensions',{'time',count(4)},'datatype','single');
        nccreate(fn_nc,'depth','Dimensions',{'layer',count(3)},'datatype','single');
        nccreate(fn_nc,'Thgt','Dimensions',{'lon',count(1),'lat',count(2), 'layer',count(3),'time',count(4)},'datatype','single');
        nccreate(fn_nc,'Tper','Dimensions',{'lon',count(1),'lat',count(2), 'layer',count(3),'time',count(4)},'datatype','single');
        nccreate(fn_nc,'Tdir','Dimensions',{'lon',count(1),'lat',count(2), 'layer',count(3),'time',count(4)},'datatype','single');
   
        hs=single(hs);
        wper=single(wper);
        wdir=single(wdir);

        xii(xii>180) = xii(xii>180)-360;

        ncwrite(fn_nc,'lon',xii);
        ncwrite(fn_nc,'lat',yii);
        ncwrite(fn_nc,'time',t);
        ncwrite(fn_nc,'depth',dep);
        ncwrite(fn_nc,'Thgt',hs);
        ncwrite(fn_nc,'Tper',wper);
        ncwrite(fn_nc,'Tdir',wdir);
    end
end