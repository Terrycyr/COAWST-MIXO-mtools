clear all;

grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v10.nc';
bay_mask0 = ncread('../Model_grid/ROMS_WFS_10river_grid_bay_mask.nc','mask_rho');
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));

out_path = strcat('./RCA/');
if(~exist(out_path))
    mkdir(out_path);
end

fn = [out_path, 'WFS_sed_2001_ini.nc'];

[r,c] = size(mask);

r=r+2;
c=c+2;

bay_mask = zeros(r,c);
bay_mask(2:end-1,2:end-1) = bay_mask0;

date_out = datenum(2001,1,1,0,0,0);
t = length(date_out);
create_rca_netcdf_init_sed(fn,gn,t);
const_initialize(fn,0.1);

CTEMP = ones(r,c)*10;

CPOPG1 = ones(r,c)*3000;
CPOPG2 = ones(r,c)*70000;
CPOPG3 = ones(r,c)*240000;

CPOPG1(bay_mask==0) = 300;
CPOPG2(bay_mask==0) = 7000;
CPOPG3(bay_mask==0) = 24000;

CPONG1 = ones(r,c)*900;
CPONG2 = ones(r,c)*18000;
CPONG3 = ones(r,c)*180000;

CPONG1(bay_mask==0) = 90;
CPONG2(bay_mask==0) = 1800;
CPONG3(bay_mask==0) = 18000;

CPOCG1 = ones(r,c)*40000;
CPOCG2 = ones(r,c)*400000;
CPOCG3 = ones(r,c)*6000000;

CPONG1(bay_mask==0) = 4000;
CPONG2(bay_mask==0) = 40000;
CPONG3(bay_mask==0) = 600000;

CPOP(:,:,1)  =CPOPG1;
CPOP(:,:,2)  =CPOPG2;
CPOP(:,:,3)  =CPOPG3;
CPON(:,:,1)  =CPONG1;
CPON(:,:,2)  =CPONG2;
CPON(:,:,3)  =CPONG3;
CPOC(:,:,1)  =CPOCG1;
CPOC(:,:,2)  =CPOCG2;
CPOC(:,:,3)  =CPOCG3;

CPOS = ones(r,c)*1300000;
CPOS(bay_mask==0) = 1300000;

%Layer 2 anaerobic
PO4T2TM1S = ones(r,c)*4400;
PO4T2TM1S(bay_mask==0) = 440;

NH4T2TM1S = ones(r,c)*18000;
NH4T2TM1S(bay_mask==0) = 1800;

NO3T2TM1S = ones(r,c)*2;

SIT2TM1S = ones(r,c)*300000;
SIT2TM1S(bay_mask==0) = 30000;

HST2TM1S = ones(r,c)*200;
CH4T2TM1S = ones(r,c)*0.5;
SO4T2TM1S = ones(r,c)*640;

BNTHSTR1S = ones(r,c)*5;

%Layer 1 aerobic
PO4T1TM1S = ones(r,c)*18000;
PO4T2TM1S(bay_mask==0) = 1800;

NH4T1TM1S = ones(r,c)*400;
NH4T1TM1S(bay_mask==0) = 40;

NO3T1TM1S = ones(r,c)*0.1;

SIT1TM1S = ones(r,c)*300000;
SIT2TM1S(bay_mask==0) = 30000;

CH4T1TM1S = ones(r,c)*0.005;
HST1TM1S = ones(r,c)*0.001;




ncwrite(fn,'CTEMP',CTEMP);
ncwrite(fn,'CPOP',CPOP);
ncwrite(fn,'CPON',CPON);
ncwrite(fn,'CPOC',CPOC);
ncwrite(fn,'CPOS',CPOS);
ncwrite(fn,'BNTHSTR1S',BNTHSTR1S);
ncwrite(fn,'PO4T2TM1S',PO4T2TM1S);
ncwrite(fn,'NH4T2TM1S',NH4T2TM1S);
ncwrite(fn,'NO3T2TM1S',NO3T2TM1S);
ncwrite(fn,'HST2TM1S',HST2TM1S);
ncwrite(fn,'SIT2TM1S',SIT2TM1S);
ncwrite(fn,'CH4T2TM1S',CH4T2TM1S);
ncwrite(fn,'SO4T2TM1S',SO4T2TM1S);
ncwrite(fn,'PO4T1TM1S',PO4T1TM1S);
ncwrite(fn,'NH4T1TM1S',NH4T1TM1S);
ncwrite(fn,'NO3T1TM1S',NO3T1TM1S);
ncwrite(fn,'HST1TM1S',HST1TM1S);
ncwrite(fn,'SIT1TM1S',SIT1TM1S);
ncwrite(fn,'CH4T1TM1S',CH4T1TM1S);

ncwrite(fn,'TIME',0);