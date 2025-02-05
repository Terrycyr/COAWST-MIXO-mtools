%Modified by terry, 2021-11-07
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

fn = [out_path,'biotest_eutr_bry.nc'];

delete(fn);

n_river = 1;
dir_flag = [1 1 1 1]; %S N W E

date_out = datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,24,0,0);
year = datevec(date_out(1));
year = year(1);
t = length(date_out);
create_rca_netcdf_bndry_eutr(fn,gn,t,n_river,dir_flag)
const_initialize(fn,0.);



%
tmp = ncread(fn,'DO_south');
DO_south = 7*ones(size(tmp));

tmp = ncread(fn,'DO_north');
DO_north = 7*ones(size(tmp));

tmp = ncread(fn,'DO_east');
DO_east = 7*ones(size(tmp));

tmp = ncread(fn,'DO_west');
DO_west = 7*ones(size(tmp));

ncwrite(fn,'bry_time',date_out-date_out(1));
ncwrite(fn,'IBCOPT',2);
ncwrite(fn,'IBCPWLOPT',1);
ncwrite(fn,'bry_time',date_out-date_out(1));

%south
ncwrite(fn,'DO_south',DO_south);

%north
ncwrite(fn,'DO_north',DO_north);

%east
ncwrite(fn,'DO_east',DO_east);

%west
ncwrite(fn,'DO_west',DO_west);
