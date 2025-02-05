clear all;
rst = './WFS_2003_hot.nc';
time = ncread(rst,'ocean_time')
%time2 = ones(size(time))*(datenum(2023,1,1)-datenum(2022,1,1))*86400;
time2 = zeros(size(time));
ncwrite(rst,'ocean_time',time2)