clear all; close all;

r_file1 = 'C:\Users\cheny\Desktop\EcoHAB\Hurricane_Ian_2022\MIXO_input\WFS_2022_2023_river_bio_mixo.nc';
r_file2 = 'WFS_2022_river_no_IAN.nc';
r_file3 = 'WFS_2022_2023_river_bio_mixo_no_IAN_nutr.nc';

delete(r_file3);
copyfile(r_file1,r_file3);

out_date1 = datenum(2022,1,1):datenum(2023,6,1);
out_date2 = datenum(2022,1,1):datenum(2023,1,1);
rflow1 = ncread(r_file1,'river_transport');
rflow2 = ncread(r_file2,'river_transport');

[r,c] = size(rflow1);
N = 21;

for i=1:length(out_date2)
    delta = abs(out_date2(i)-out_date1);
    out_i(i) = find(delta==min(delta));
end

Nutri_adjust = ones(r,21,c);
for i=1:N
    Nutri_adjust(:,i,out_i) = rflow2./rflow1(:,out_i);
end

Nutri_adjust(isnan(Nutri_adjust)) = 1;

river_LDOC = ncread(r_file3,'river_LDOC');
river_LDON = ncread(r_file3,'river_LDON');
river_LDOP = ncread(r_file3,'river_LDOP');
river_RDOC = ncread(r_file3,'river_RDOC');
river_RDON = ncread(r_file3,'river_RDON');
river_RDOP = ncread(r_file3,'river_RDOP');
river_NH4T = ncread(r_file3,'river_NH4T');
river_NO23 = ncread(r_file3,'river_NO23');
river_PO4T = ncread(r_file3,'river_PO4T');
river_LPOC = ncread(r_file3,'river_LPOC');
river_LPON = ncread(r_file3,'river_LPON');
river_LPOP = ncread(r_file3,'river_LPOP');
river_RPOC = ncread(r_file3,'river_RPOC');
river_RPON = ncread(r_file3,'river_RPON');
river_RPOP = ncread(r_file3,'river_RPOP');
river_SIT = ncread(r_file3,'river_SIT');

ncwrite(r_file3,'river_LDOC',river_LDOC.*Nutri_adjust);
ncwrite(r_file3,'river_LDON',river_LDON.*Nutri_adjust);
ncwrite(r_file3,'river_LDOP',river_LDOP.*Nutri_adjust);
ncwrite(r_file3,'river_RDOC',river_RDOC.*Nutri_adjust);
ncwrite(r_file3,'river_RDON',river_RDON.*Nutri_adjust);
ncwrite(r_file3,'river_RDOP',river_RDOP.*Nutri_adjust);
ncwrite(r_file3,'river_NH4T',river_NH4T.*Nutri_adjust);
ncwrite(r_file3,'river_NO23',river_NO23.*Nutri_adjust);
ncwrite(r_file3,'river_PO4T',river_PO4T.*Nutri_adjust);
ncwrite(r_file3,'river_LPOC',river_LPOC.*Nutri_adjust);
ncwrite(r_file3,'river_LPON',river_LPON.*Nutri_adjust);
ncwrite(r_file3,'river_LPOP',river_LPOP.*Nutri_adjust);
ncwrite(r_file3,'river_RPOC',river_RPOC.*Nutri_adjust);
ncwrite(r_file3,'river_RPON',river_RPON.*Nutri_adjust);
ncwrite(r_file3,'river_RPOP',river_RPOP.*Nutri_adjust);
ncwrite(r_file3,'river_SIT',river_SIT.*Nutri_adjust);





