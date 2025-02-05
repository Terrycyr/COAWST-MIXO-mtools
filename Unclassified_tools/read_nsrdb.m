clear all; close all;
fname = '1005027_27.05_-82.74_2002.csv';
[dat,txt,raw] = xlsread(fname);
year = dat(3:end,1);
mon = dat(3:end,2);
day = dat(3:end,3);
HH = dat(3:end,4);
MM = dat(3:end,5);
%dni = dat(3:end,8);
ghi = dat(3:end,8);
wind = dat(3:end,17);

time = datenum(year,mon,day,HH,MM,0);
date_out = datenum(2002,1,1,0,0,0):1/24:datenum(2002,12,31,24,0,0);

ghi_out = interp1(time,ghi,date_out,'linear','extrap');
wind_out = interp1(time,wind,date_out,'linear','extrap');

ITOTSF = ghi_out/0.484583; %W/m2 -> ly/day
WIND = wind_out;

save('NSRDB.mat','ITOTSF','WIND');

