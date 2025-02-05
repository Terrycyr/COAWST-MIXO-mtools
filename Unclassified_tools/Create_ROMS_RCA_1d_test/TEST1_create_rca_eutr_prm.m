clear all;

grd_name =  './bio1dtest_grd.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));

out_path = strcat('./RCA/');
if(~exist(out_path))
    mkdir(out_path);
end

fn = [out_path, 'biotest_eutr_prm.nc'];

[r,c] = size(mask);

date_out = datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,24,0,0);
date_out2 = datenum(2001,1,1,0,0,0):1/24:datenum(2001,12,31,24,0,0);
t = length(date_out);
t2 = length(date_out2);

time = date_out-date_out(1);
time2 = date_out2-date_out2(1);
create_rca_netcdf_prm_eutr(fn,gn,t,t,t2,t,t,t)
const_initialize(fn,0);

VSNET1 = ones(r+2,c+2).*0.0;
VSNET2 = ones(r+2,c+2).*0.0;
VSNET3 = ones(r+2,c+2).*0.0;
VSNET4 = ones(r+2,c+2).*0.0;

ITOTSF = ones(366,1).*100.0;
F = ones(366,1).*0.5;

ncwrite(fn,'WIND_time',time2);
ncwrite(fn,'ITOTSF_time',time);
ncwrite(fn,'ITOTSF',ITOTSF);
ncwrite(fn,'F_time',time);
ncwrite(fn,'F',F);
ncwrite(fn,'KETVF_time',time);
ncwrite(fn,'TSS_time',time);
ncwrite(fn,'Turb_time',time);
ncwrite(fn,'VSNET1',VSNET1);
ncwrite(fn,'VSNET2',VSNET2);
ncwrite(fn,'VSNET3',VSNET3);
ncwrite(fn,'VSNET4',VSNET4);
