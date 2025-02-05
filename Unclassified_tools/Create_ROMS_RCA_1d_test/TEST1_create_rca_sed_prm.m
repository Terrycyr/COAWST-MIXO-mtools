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

fn = [out_path, 'biotest_sed_prm.nc'];

create_rca_netcdf_prm_sed(fn,gn)
zero_initialize(fn);

[r,c] = size(mask);


HSED = ones(r+2,c+2)*15;
VSED = ones(r+2,c+2)*0.25;
VDMIX = ones(r,c)*0.001;
VPMIX = ones(r,c)*0.00012;

ncwrite(fn,'HSED',HSED);
ncwrite(fn,'VSED',VSED);
ncwrite(fn,'VPMIX',VPMIX);
ncwrite(fn,'VDMIX',VDMIX);