function [out_ij] = get_station_ij(location)
%get_station_ij Summary of this function goes here
%   Detailed explanation goes here
out_ij = [];
fid = fopen("ps_input",'rt+');
flag = 0;
k = 0;

while(~feof(fid))
    dat = fgetl(fid);
    if(contains(dat,'Location'))
        flag = 1;
        continue;
    end

    if(contains(dat,'-99'))
        if(k<length(location))
            no_data = find(out_ij(:,1)==0);
            no_data_st = location(no_data);
        end
        disp('No data found for station:');
        for i = 1:length(no_data)
            disp(no_data_st{i});
        end
        break;
    end

    if(flag==1)

        tmp = str2num(dat(1:16));
        st_name = strtrim(dat(17:end));
        for i=1:length(location)
            if(contains(st_name,location{i},'IgnoreCase',true))
                k = k+1;
                out_ij(i,1) = tmp(2);
                out_ij(i,2) = tmp(3);
            end
        end
    end

    if(k==length(location))
        break;
    end
end
end