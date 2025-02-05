clear all;

grd_name =  './bio1dtest_grd.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N = 21;

%DISSOLVED & PARTICULATE
OCDP_r = 0.84;
ONDP_r = 0.60;
OPDP_r = 0.22;

%LABILE & REFRACTORY
OCLR_r = 0.35;
ONLR_r = 0.35;
OPLR_r = 0.35;

%C/N,C/P,C/Si
CRBP1 = 40.;
CRBN1 = 5.67;
CRBS1 = 3.;
CRBP2 = 32.;
CRBN2 = 6.3;
CRBS2 = 1e21;
CRBP3 = 62.;
CRBN3 = 4.67;
CRBS3 = 1e21;

out_path = strcat('./RCA/');
if(~exist(out_path,"dir"))
    mkdir(out_path);
end

fn = [out_path, 'biotest_eutr_ini.nc'];

delete(fn);

date_out = datenum(2001,1,1,0,0,0);
year = datevec(date_out(1));
year = year(1);
t = length(date_out);
create_rca_netcdf_init_eutr(fn,gn,t)
const_initialize(fn,0);

%
ncwrite(fn,'TIME',0);
ncwrite(fn,'SAL',0*repmat(ones(size(lon)),1,1,21));
ncwrite(fn,'DO',7*repmat(ones(size(lon)),1,1,21));