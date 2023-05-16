function river_out = get_nwis_timeseries(r,out_date,varname)
%get_nwis_timeseries 此处显示有关此函数的摘要
%   此处显示详细说明

    for n_r = 1:size(r,2)
        river_out_m = zeros(1,12);
        river_m_count = zeros(1,12);
        river_date0 = [];
        river_out0 = [];
        x = strfind(r(n_r).descriptions,'Calendar');
        tmp = find(cellfun(@(x)~isempty(x),x,'UniformOutput',true));
        pos_year = tmp(1);pos_mon = tmp(2);pos_day = tmp(3);
        x = strfind(r(n_r).descriptions,varname);
        tmp = find(cellfun(@(x)~isempty(x),x,'UniformOutput',true));
        try
            pos_out = tmp(1);
        catch
            river_out(n_r,:) = NaN*ones(1,length(out_date));
            continue;
        end          

        for i=1:size(r(n_r).data,1)
            river_date0(i) = datenum(r(n_r).data(i,pos_year),...
            r(n_r).data(i,pos_mon),r(n_r).data(i,pos_day))+5/24; %convert to UTC
            river_out0(i) = r(n_r).data(i,pos_out);
            if(~isnan(river_out0(i)))     
                river_out_m(r(n_r).data(i,pos_mon)) = river_out_m(r(n_r).data(i,pos_mon))...
                    +river_out0(i);
                river_m_count(r(n_r).data(i,pos_mon)) = river_m_count(r(n_r).data(i,pos_mon))+1;
            end
        end

        for mon=1:12
            river_out_m(mon) = river_out_m(mon)/river_m_count(mon);
        end
        river_out(n_r,:) = interp1(river_date0(~isnan(river_out0)),river_out0(~isnan(river_out0)),out_date);
        river_out_cycle = zeros(1,length(out_date));
        for t = 1:length(out_date)
            mon = str2double(datestr(out_date(t),'mm'));
            river_out_cycle(t) = river_out_m(mon);
        end
        river_out(n_r,isnan(river_out(n_r,:))) = river_out_cycle(isnan(river_out(n_r,:)));
    end

    
    

    

end

