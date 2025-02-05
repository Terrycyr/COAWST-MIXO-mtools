clear all;
hot = './WFS_2001_hot.nc';
dye01 = ncread('../NC_file_generation/WFS_2001_ini.nc','dye_01');
ncwrite(hot, 'dye_01',dye01);