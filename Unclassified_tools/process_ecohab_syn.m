clear all;

[~,~,raw] = xlsread('ECOHAB-Karenia Data -4 blooms from Access-WORKING COMBO FILE (ALL STATIONS) 3-31-23.xlsx',2);

k=0;

for i=4:764
    k=k+1;
    dat_date(k) = datenum(raw{i,8},'mm/dd/yyyy');
    dat_lat(k) = raw{i,5};
    dat_lon(k) = raw{i,6};
    dat_dep(k) = raw{i,7};
    dat_kb(k) = raw{i,26};
    dat_syn(k) = raw{i,124};
    dat_sal(k) = raw{i,43};
    dat_temp(k) = raw{i,45};
    %uM
    dat_dpo4(k) = raw{i,64};
    dat_dnh4(k) = raw{i,65};
    dat_dnox(k) = raw{i,65};
    dat_tdp(k) = raw{i,67};
    dat_tdn(k) = raw{i,68};
    dat_don(k) = raw{i,82};
    dat_dop(k) = raw{i,83};
    dat_tn(k) = raw{i,84};
    dat_tp(k) = raw{i,85};
    dat_sio2(k) = raw{i,69};
end

save('OLD_ECOHAB_dat.mat',"dat_date","dat_lat","dat_lon","dat_dep"...
    ,"dat_kb","dat_syn","dat_sal","dat_temp","dat_dpo4","dat_dnh4","dat_dnox","dat_tdp"...
    ,"dat_tdn","dat_tn","dat_tp","dat_sio2","-v7.3");