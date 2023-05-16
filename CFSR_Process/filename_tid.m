function [cfsr_fname,tid] = filename_tid(cfsr_dir,time,varname)
%filename_tid Summary of this function goes here
%   Detailed explanation goes here
    
    if(time>datenum(2011,4,1))
        tmp = datevec(time);
        hour = tmp(4);
        if(hour==0)
            cfsr_fname0=dir(strcat(cfsr_dir,'*',datestr(time-1,'yyyymmdd'),'*.nc'));
           
        else
            cfsr_fname0=dir(strcat(cfsr_dir,'*',datestr(time,'yyyymmdd'),'*.nc'));
        end
        cfsr_fname = strcat(cfsr_dir,cfsr_fname0.name);
        ncid = netcdf.open(cfsr_fname,'NOWRITE');
        varid = netcdf.inqVarID(ncid,varname);
        [~,~,dimid,~] = netcdf.inqVar(ncid,varid);
        [time_name,~] = netcdf.inqDim(ncid,dimid(end));
        cfsr_hour = ncread(cfsr_fname,time_name);
        cfsr_hour(cfsr_hour==24) = 0; 
    else
        file_list = dir(strcat(cfsr_dir,'*','-','*.nc'));
        %date list
        for i = 1:length(file_list)
            fname{i} = file_list(i).name;
            bar_loc = strfind(fname{i},'-');
            tmp = fname{i}([-8:-1]+bar_loc);
            date_list(i,1) = datenum(str2double(tmp(1:4)),str2double(tmp(5:6)),str2double(tmp(7:8)));
            tmp = fname{i}([1:8]+bar_loc);
            date_list(i,2) = datenum(str2double(tmp(1:4)),str2double(tmp(5:6)),str2double(tmp(7:8)),24,0,0);
            f_flag(i) = (time>date_list(i,1))&&(time<=date_list(i,2));
        end
        f_pos = find(f_flag);
        hour = (time-date_list(f_pos,1))*24;
        if(hour==0)
            cfsr_fname = strcat(cfsr_dir,fname{f_pos-1});
            ncid = netcdf.open(cfsr_fname,'NOWRITE');
            varid = netcdf.inqVarID(ncid,varname);
            [~,~,dimid,~] = netcdf.inqVar(ncid,varid);
            [time_name,~] = netcdf.inqDim(ncid,dimid(end));
            cfsr_hour = ncread(cfsr_fname,time_name);
            cfsr_hour(end) = 0;
        else
            cfsr_fname = strcat(cfsr_dir,fname{f_pos});
            ncid = netcdf.open(cfsr_fname,'NOWRITE');
            varid = netcdf.inqVarID(ncid,varname);
            [~,~,dimid,~] = netcdf.inqVar(ncid,varid);
            [time_name,~] = netcdf.inqDim(ncid,dimid(end));
            cfsr_hour = ncread(cfsr_fname,time_name);
        end
        
    end
    netcdf.close(ncid);
    [~,tid] = min(abs(hour-cfsr_hour));
    
end

