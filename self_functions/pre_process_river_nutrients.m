function [nutri_dat,nutri_mean] = pre_process_river_nutrients(fn,varname,n_river,on_flag)
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
for i=1:n_river
    tmp = xlsread(fn,i);
    k0=0;
    for j=1:length(varname)
        if(strcmp(varname{j},'OP'))
            continue;
        end
        k0=k0+1;
        dat_tmp = tmp(:,2*k0);
        t_tmp = xlstime2date(tmp(:,2*k0-1));
        t_tmp(isnan(dat_tmp)) = [];
        dat_tmp(isnan(dat_tmp)) = [];
        if(length(unique(t_tmp))<12)
            t_tmp = t_tmp+rand(length(t_tmp),1);
        end
        dat_raw{i} = dat_tmp;

        if(~isempty(dat_tmp))
            time_min(i) = min(t_tmp);
            time_max(i) = max(t_tmp);
            time_unique = unique(t_tmp);

            for k=1:length(time_unique)
                t(k) = time_unique(k);
                value(k) = mean(dat_tmp(t_tmp==time_unique(k)));
            end

            nutri_dat{i,j}(:,1) = t;
            nutri_dat{i,j}(:,2) = value;
            nutri_dat{i,j} = sortrows(nutri_dat{i,j},1);

            time2 = time_min(i):time_max(i);
            
            tmp2 = interp1(nutri_dat{i,j}(:,1),nutri_dat{i,j}(:,2),time2);
            timev = datevec(time2);
            time_flag = timev(:,2)+timev(:,3)/100;
            time_unique2 = unique(time_flag);

            for k=1:length(time_unique2)
                t2(k) = time_unique2(k);
                value2(k) = mean(tmp2(time_flag==time_unique2(k)));
            end

            nutri_mean{i,j}(:,1) = t2;
            nutri_mean{i,j}(:,2) = smooth(value2,7);

            if(strcmp(varname{j},'TN')&&on_flag(i)==1)
                nh4_tmp = interp1(nutri_dat{i,1}(:,1),nutri_dat{i,1}(:,2),nutri_dat{i,j}(:,1));
                no23_tmp = interp1(nutri_dat{i,2}(:,1),nutri_dat{i,2}(:,2),nutri_dat{i,j}(:,1));
                on_tmp = nutri_dat{i,j}(:,2)-nh4_tmp-no23_tmp;
                t_tmp2(:,1) = nutri_dat{i,j}(:,1);
                t_tmp2(on_tmp<0|isnan(on_tmp)) = [];
                on_tmp(on_tmp<0|isnan(on_tmp)) = [];
                nutri_dat{i,3}(:,1) = t_tmp2;
                nutri_dat{i,3}(:,2) = on_tmp;
            end

            clear t value t2 value2 t_tmp t_tmp2
        end
    end
    po4_tmp = interp1(nutri_dat{i,5}(:,1),nutri_dat{i,5}(:,2),nutri_dat{i,6}(:,1));
    op_tmp = nutri_dat{i,6}(:,2)-po4_tmp;
    t_tmp(:,1) = nutri_dat{i,6}(:,1);
    t_tmp(op_tmp<0|isnan(op_tmp)) = [];
    op_tmp(op_tmp<0|isnan(op_tmp)) = [];
    nutri_dat{i,9}(:,1) = t_tmp;
    nutri_dat{i,9}(:,2) = op_tmp;
end

end
