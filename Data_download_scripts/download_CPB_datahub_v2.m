clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%station_list
% station_name={'CB3.3E','CB3.3C','CB3.3W','CB4.4','CB5.4','CB8.1','LE2.3','LE5.5A','LE5.5B'};
% station_name={'CB2.2','CB3.1','CB3.2','CB3.3W','CB3.3C','CB3.3E',...
%     'CB4.1W','CB4.1C','CB4.1E','CB4.2W','CB4.2C','CB4.2E','CB4.3W',...
%     'CB4.3C','CB4.3E','CB4.4','CB5.1W','CB5.1','CB5.2','CB5.3','CB5.4',...
%     'CB5.5','CB7.1N','CB7.1','CB6.1','CB7.1S','CB6.2','CB6.3','CB7.2',...
%     'CB7.2E','CB6.4','CB7.3E','CB7.3','CB7.4N','CB7.4','CB8.1','CB8.1E',...
%     'LE2.3','LE3.6','LE4.3','LE5.5A','LE5.5B'};

if(1)
    [~,~,station_name]=xlsread('all_station.xlsx')
    station_name=station_name';
else
    station_name={'CB3.2'};
end

%variable_list
vname={'SALINITY','WTEMP','DO'};

%Output directory
out_dir = './';

for iyear=2003

    fid=fopen(['CPB_monitoring_stations_lonlat_',num2str(iyear),'.dat'],'w');
    fid2=fopen(['CPB_monitoring_stations_name_',num2str(iyear),'.dat'],'w');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %no change below
    Start_Date=['1-1-',num2str(iyear)];
    End_Date=['1-1-',num2str(iyear+1)];

    ii=0;

    for istation=1:size(station_name,2)
        for ivar=1:size(vname,2)

            disp(sprintf('%-10s%10s%10s',vname{ivar},station_name{istation},num2str(iyear)))
            Data_Stream_Value ='0,1';

            url='https://data.chesapeakebay.net/api.json/WaterQuality/Programs';
            a1=webread(url);
            a2=struct2cell(a1);
            a3=cell2mat(a2(1,:));
            a4=sprintf('%d,',a3);
            Program_Id=a4(1:end-1);

            url='https://data.chesapeakebay.net/api.json/WaterQuality/Projects';
            a1=webread(url);
            a2=struct2cell(a1);
            a3=cell2mat(a2(1,:));
            a4=sprintf('%d,',a3);
            Project_Id=a4(1:end-1);

            url='https://data.chesapeakebay.net/api.json/Station';
            a1=webread(url);
            a2=struct2cell(a1);
            a3=a2(2,:);
            a4=find(strcmp(a3,station_name{istation}));
            a5=a2{1,a4};
            Attribute_Id=num2str(a5);

            url='https://data.chesapeakebay.net/api.json/Substances';
            a1=webread(url);
            a2=struct2cell(a1);
            a3=a2(2,:);
            a4=find(strcmp(a3,vname{ivar}));
            a5=a2{1,a4};
            Substance_Id=num2str(a5);

            namespace='https://data.chesapeakebay.net/api.json/WaterQuality/WaterQuality';

            url=sprintf('%s/',namespace,Start_Date,End_Date,Data_Stream_Value,Program_Id,Project_Id,'Station',Attribute_Id,Substance_Id);

            dat0=webread(url);
            try
                dat=struct2cell(dat0);
            catch
                disp('------>>>>> no data <<<<<------')
                continue
            end

            sample_date=dat(9,:);
            sample_time=dat(10,:);
            obs_value=cell2mat(dat(20,:))';
            unit=dat(21,:)';
            total_depth=cell2mat(dat(11,:))';
            obs_depth=cell2mat(dat(14,:))';
            obs_layer=dat(15,:)';
            obs_lon=cell2mat(dat(29,:))';
            obs_lat=cell2mat(dat(28,:))';

            nt=size(sample_time,2);
            for i=1:nt
                t1=sample_date{1,i};
                t2=sample_time{1,i};
                t3=split(t2,':');
                t4=str2num(t3{1})/24+str2num(t3{2})/1440+str2num(t3{3})/86400;
                obs_time(i,1)=datenum(sample_date(1,i),'yyyy-mm-ddTHH:MM:SS')+t4;
            end

            obs_Times=datevec(obs_time);
            output=[out_dir,vname{ivar},'_',station_name{istation},'_',num2str(iyear),'.mat'];
            save(output,"obs_time","obs_Times","total_depth","obs_depth","obs_layer","obs_lon","obs_lat","obs_value","unit")

            if(ivar==1)
                fprintf(fid,'%12.4f %12.4f \n',obs_lon(1),obs_lat(1));
                fprintf(fid2,'%20s \n',station_name{istation});

                ii=ii+1;
                name{ii,:}=station_name{istation};
            end
            clear obs_time obs_Times total_depth obs_depth obs_layer obs_lon obs_lat obs_value unit
        end
    end
    fclose(fid);
    fclose(fid2);

    save(['CPB_monitoring_stations_name_',num2str(iyear),'.mat'],'name')
    clear name
end

tada