clear all;

[~,~,raw] = xlsread('./SYNCOC_FINAL.xlsx');

k = 0;
for i = 2:size(raw,1)
    k=k+1;
    s_lon(k) = raw{i,5};
    s_lat(k) = raw{i,4};
    s_year(k) = raw{i,2};
    s_name{k} = raw{i,3};
    s_syn(k) = raw{i,6}*1000;
end

save('ECOHAB_syn_dat.mat',"s_year","s_lat","s_lon","s_name","s_syn","-v7.3");