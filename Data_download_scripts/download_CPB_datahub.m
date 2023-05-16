clear;clc;tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%station_list
station_name={'CB3.3E','CB3.3C','CB3.3W','CB4.4','CB5.4','CB8.1','LE2.3','LE5.5A','LE5.5B'};
%variable_list
vname={'SALINITY','WTEMP'};
year=[2003,2004];

for iyear=year
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %no change below
    Start_Date=['1-1-',num2str(iyear)];
    End_Date=['1-1-',num2str(iyear+1)];

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
            value=cell2mat(dat(20,:))';
            unit=dat(21,:)';
            total_depth=cell2mat(dat(11,:))';
            obs_depth=cell2mat(dat(14,:))';
            obs_lon=cell2mat(dat(29,:))';
            obs_lat=cell2mat(dat(28,:))';


            nt=size(sample_time,2);
            for i=1:nt
                t1=sample_date{1,i};
                t2=sample_time{1,i};
                t3=split(t2,':');
                t4=str2num(t3{1})/24+str2num(t3{2})/1440+str2num(t3{3})/86400;
                time(i,1)=datenum(sample_date(1,i),'yyyy-mm-ddTHH:MM:SS')+t4;
            end

            Times=datevec(time);
            output=[vname{ivar},'_',station_name{istation},'_',num2str(iyear),'.mat'];
            save(output,"time","Times","total_depth","obs_depth","obs_lon","obs_lat","value","unit")

        end
    end
end
toc
disp("TADA")