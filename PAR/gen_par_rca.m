clear all;
load erdMGpar01day_c46c_e27d_dbdf.mat

out_date1 = datenum(2001,1,1,0,0,0):3/24:datenum(2001,12,31,24,0,0);

time = erdMGpar01day.time;
tvec = datevec(time/3600/24+datenum(1970,1,1));
par = mean(mean(double(erdMGpar01day.par),4),3);
for i = 1:366
    tmp = datevec(i-2+datenum(1970,1,1));
    mon = tmp(2);
    day = tmp(3);
    pos = find(tvec(:,2)==mon&tvec(:,3)==day);
    date(i) = mean(time(pos)/3600/24+datenum(1970,1,1));
    par_m(i) = mean(par(pos),'omitnan')*1e6/24/3600; %Einsteins m-2 d-1 -> µmoles m-2 s-1
end

ITOTSF = par_m;

save('PAR_rca.mat','ITOTSF');