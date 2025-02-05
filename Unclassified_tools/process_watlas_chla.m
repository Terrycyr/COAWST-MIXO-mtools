clear all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

[~,~,raw] = xlsread('./WaterAtlas_CHLA2.xlsx');

k = 0;
for i = 2:size(raw,1)
    if(sum(isnan(raw{i,18}))&&~isnan(raw{i,16}))
        if(~isempty(strfind(raw{i,14},'Chlorophyll a')))
            k=k+1;
            s_lon(k) = raw{i,8};
            s_lat(k) = raw{i,7};
            s_date(k) = xlstime2date(raw{i,10});
            if(strcmp(raw{i,12},'m'))
                s_dep(k) = raw{i,11};
            elseif(strcmp(raw{i,12},'ft'))
                s_dep(k) = raw{i,11}*0.3048;
            end
            s_chla(k) = raw{i,16};
        end
    end
end

save('WA_chla_dat.mat',"s_date","s_lat","s_lon","s_dep","s_chla","-v7.3");