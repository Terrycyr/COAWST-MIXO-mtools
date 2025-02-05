clear all;

rst = './WFS_2023_no_ian_hot.nc';
time = ncread(rst,'ocean_time')
time(1) = (datenum(2023,1,1)-datenum(2022,1,1))*86400;
ncwrite(rst,'ocean_time',time);