function [kb_lon,kb_lat, kb_dep, kb_dat, kb_dep_all, kb_dat_all] = get_FWC_kb(t_min,t_max,mat_data,method)
%get_FWC_kb Summary of this function goes here
%   Detailed explanation goes here

    load(mat_data);

    pos = s_date>=t_min&s_date<=t_max;

    lon_tmp = s_lon(pos);
    lat_tmp = s_lat(pos);
    dep_tmp = s_dep(pos);
    kb_tmp = s_kb(pos);

    if(strcmp(method,'mean'))
        [C,ia,ic] = unique([lon_tmp',lat_tmp'],'rows','stable');
        for i=1:length(ia)
            kb_lon(i) = C(i,1);
            kb_lat(i) = C(i,2);
            kb_dep(i) = mean(dep_tmp(ic==i),'omitnan');
            kb_dep_all{i} = dep_tmp(ic==i);
            kb_dat(i) = mean(kb_tmp(ic==i),'omitnan');
            kb_dat_all{i} = kb_tmp(ic==i);
        end
    end

end