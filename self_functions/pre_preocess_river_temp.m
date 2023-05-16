function [temp_dat,temp_mean] = pre_preocess_river_temp(fn,n_river)

for i=1:n_river
    tmp = xlsread(fn,i);
    tmp(:,1) = datenum(1900,1,1)+tmp(:,1)+5/24;
    temp_raw{i} = tmp;
    time_min(i) = min(tmp(:,1));
    time_max(i) = max(tmp(:,1));
    time_unique = unique(tmp(:,1));

    for j=1:length(time_unique)
        t(j) = time_unique(j);
        value(j) = mean(tmp(tmp(:,1)==time_unique(j),2));
    end
    temp_dat{i}(:,1) = t;
    temp_dat{i}(:,2) = value;
    temp_dat{i} = sortrows(temp_dat{i},1);

    time2 = time_min(i):time_max(i);
    tmp2 = interp1(temp_dat{i}(:,1),temp_dat{i}(:,2),time2);
    timev = datevec(time2);
    time_flag = timev(:,2)+timev(:,3)/100;
    time_unique2 = unique(time_flag);

    for j=1:length(time_unique2)
        t2(j) = time_unique2(j);
        value2(j) = mean(tmp2(time_flag==time_unique2(j)));
    end

    temp_mean{i}(:,1) = t2;
    temp_mean{i}(:,2) = smooth(value2,7);

    clear t value t2 value2
end

end

