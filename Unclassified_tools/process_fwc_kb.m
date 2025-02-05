clear all;

[~,~,raw] = xlsread('./Historic_Harmful_Algal_Bloom_Events_2000_-_2006.csv');

k = 0;
for i = 2:size(raw,1)
    if(strcmp(raw{i,11},'Karenia brevis'))
        k=k+1;
        s_lon(k) = raw{i,10};
        s_lat(k) = raw{i,9};
        s_date(k) = datenum(raw(i,4),'yyyy/mm/dd HH:MM:SS');
        s_dep(k) = raw{i,7};
        s_kb(k) = raw{i,12};
    end
end

save('FWC_kb_dat_0006.mat',"s_date","s_lat","s_lon","s_dep","s_kb","-v7.3");