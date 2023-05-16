clear all;

grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));

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

fn = [out_path, 'WFS_eutr_2001_ini.nc'];

delete(fn);

date_out = datenum(2001,1,1,0,0,0);
year = datevec(date_out(1));
year = year(1);
t = length(date_out);
create_rca_netcdf_init_eutr(fn,gn,t)
const_initialize(fn,0);

load(strcat('../Water_Atlas/','bay_ini_',num2str(year),'.mat'));
load(strcat('../GOM_preprocessing/','ocean_nutrients_ini_',num2str(year),'.mat'));

%
SAL = sal;
DO = do;
NO23 = no3/10;
NH4T = nh4;
PO4T = po4;
SIT = si;
LDON = zeros(size(si));
RDON = zeros(size(si));
LPON = zeros(size(si));
RPON = zeros(size(si));
LDOC = zeros(size(si));
RDOC = zeros(size(si));
LPOC = zeros(size(si));
RPOC = zeros(size(si));

%flag = repmat(h<40&mask==1,1,1,size(si,3));
%flag(410:end,:,:) = 0;
flag = 1;
PHYT1 = ones(size(si))*1e-6.*flag;
PHYT2 = ones(size(si))*1e-7.*flag;
PHYT3 = ones(size(si))*1e-6.*flag;

for i=1:gn.N
    tmp = SAL(:,:,i);
    tmp2 = sal_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    SAL(:,:,i) = tmp;

    tmp = NO23(:,:,i);
    tmp2 = no23_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    NO23(:,:,i) = tmp;
    
    tmp = PO4T(:,:,i);
    tmp2 = po4_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    PO4T(:,:,i) = tmp;
    
    tmp = SIT(:,:,i);
    tmp2 = sit_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    SIT(:,:,i) = tmp;
    
    tmp = NH4T(:,:,i);
    tmp2 = nh4_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    NH4T(:,:,i) = tmp;

    tmp = LDON(:,:,i);
    tmp2 = don_bay(:,:,i);

    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*ONLR_r;
    LDON(:,:,i) = tmp;

    tmp = RDON(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONLR_r);
    RDON(:,:,i) = tmp;

    tmp = LPON(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)/ONDP_r*ONLR_r;
    LPON(:,:,i) = tmp;

    tmp = RPON(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)/ONDP_r*(1-ONLR_r);
    RPON(:,:,i) = tmp;

    tmp = LDOC(:,:,i);
    tmp2 = oc_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*ONDP_r*ONLR_r;
    LDOC(:,:,i) = tmp;

    tmp = RDOC(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*ONDP_r*(1-ONLR_r);
    RDOC(:,:,i) = tmp;

    tmp = LPOC(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)*ONLR_r;
    LPOC(:,:,i) = tmp;

    tmp = RPOC(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)*(1-ONLR_r);
    RPOC(:,:,i) = tmp;

    tmp = PHYT1(:,:,i);
    tmp2 = chla_bay(:,:,i)*50*0.5;
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    PHYT1(:,:,i) = tmp;
end

SIT = SIT+PHYT1*(1/CRBS1)+PHYT2*(1/CRBS2)+PHYT3*(1/CRBS3);
PO4T = PO4T+PHYT1*(1/CRBP1)+PHYT2*(1/CRBP2)+PHYT3*(1/CRBP3);
NH4T = NH4T+PHYT1*(1/CRBN1)+PHYT2*(1/CRBN2)+PHYT3*(1/CRBN3);

[a,b,c,d] = size(SAL);
%
ncwrite(fn,'TIME',0);
ncwrite(fn,'SAL',SAL,[2 2 1 1]);
ncwrite(fn,'SIT',SIT,[2 2 1 1]);
ncwrite(fn,'DO',DO,[2 2 1 1]);
ncwrite(fn,'NO23',NO23,[2 2 1 1]);
ncwrite(fn,'PO4T',PO4T,[2 2 1 1]);
ncwrite(fn,'NH4T',NH4T,[2 2 1 1]);
ncwrite(fn,'PHYT1',PHYT1,[2 2 1 1]);
ncwrite(fn,'PHYT2',PHYT2,[2 2 1 1]);
ncwrite(fn,'PHYT3',PHYT3,[2 2 1 1]);
ncwrite(fn,'LDON',LDON,[2 2 1 1]);
ncwrite(fn,'RDON',RDON,[2 2 1 1]);
ncwrite(fn,'LPON',LPON,[2 2 1 1]);
ncwrite(fn,'RPON',RPON,[2 2 1 1]);
ncwrite(fn,'LDOC',LDOC,[2 2 1 1]);
ncwrite(fn,'RDOC',RDOC,[2 2 1 1]);
ncwrite(fn,'LPOC',LPOC,[2 2 1 1]);
ncwrite(fn,'RPOC',RPOC,[2 2 1 1]);
