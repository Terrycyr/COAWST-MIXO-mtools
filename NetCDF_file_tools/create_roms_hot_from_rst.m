clear all;
rst = './WFS_2021_hot.nc';
time = ncread(rst,'ocean_time')
bay_mask = ncread('../Model_grid/ROMS_WFS_10river_grid_bay_mask.nc','mask_rho');

salt = ncread('../NC_file_generation/WFS_2021_ini.nc','salt');
salt2 = ncread(rst,'salt');

dye01 = ncread('../NC_file_generation/WFS_2021_ini.nc','dye_01');

temp = ncread('../NC_file_generation/WFS_2021_ini.nc','temp');

for i = 1:size(salt2,3)
    for j=1:size(salt2,4)
        tmp = salt2(:,:,i,j);
        tmp2 = salt(:,:,i);
        tmp(bay_mask==1) = tmp2(bay_mask==1);
        salt2_new(:,:,i,j) = tmp;
    end
end

ncwrite(rst,'temp',temp);
ncwrite(rst,'salt',salt2_new);
ncwrite(rst,'dye_01',dye01);

ncwrite(rst,'ocean_time',zeros(size(time)))