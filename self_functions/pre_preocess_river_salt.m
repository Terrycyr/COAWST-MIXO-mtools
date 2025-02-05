function [salt_dat,salt_mean] = pre_preocess_river_salt(fn,n_river)
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
for i=1:n_river
    tmp = xlsread(fn,i);
    if(~isempty(~isnan(tmp)))
        tmp(:,1) = xlstime2date(tmp(:,1));
        salt_raw{i} = tmp;
        time_min(i) = min(tmp(:,1));
        time_max(i) = max(tmp(:,1));
        time_unique = unique(tmp(:,1));

        for j=1:length(time_unique)
            t(j) = time_unique(j);
            value(j) = mean(tmp(tmp(:,1)==time_unique(j),2));
        end
        salt_dat{i}(:,1) = t;
        salt_dat{i}(:,2) = value;
        salt_dat{i} = sortrows(salt_dat{i},1);

        time2 = time_min(i):time_max(i);
        tmp2 = interp1(salt_dat{i}(:,1),salt_dat{i}(:,2),time2);
        timev = datevec(time2);
        time_flag = timev(:,2)+timev(:,3)/100;
        time_unique2 = unique(time_flag);

        for j=1:length(time_unique2)
            t2(j) = time_unique2(j);
            value2(j) = mean(tmp2(time_flag==time_unique2(j)));
        end

        salt_mean{i}(:,1) = t2;
        salt_mean{i}(:,2) = smooth(value2,7);

        clear t value t2 value2
    else
        salt_dat{i} = [];
        salt_mean{i} = [];
    end
end

end

