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

fn = [out_path, 'biotest_sed_ini.nc'];

[r,c] = size(mask);

r=r+2;
c=c+2;

date_out = datenum(2001,1,1,0,0,0);
t = length(date_out);
create_rca_netcdf_init_sed(fn,gn,t);
const_initialize(fn,0);



%DISSOLVED & PARTICULATE, Riverine
OCDP_r = 0.84;
ONDP_r = 0.60;
OPDP_r = 0.22;

%LABILE & REFRACTORY, Riverine
OCLR_r = 0.35;
ONLR_r = 0.35;
OPLR_r = 0.35;

%DISSOLVED & PARTICULATE, Offshore
OCDP_r_off = 0.84;
ONDP_r_off = 0.60;
OPDP_r_off = 0.22;

%LABILE & REFRACTORY, Offshore
OCLR_r_off = 0.35;
ONLR_r_off = 0.35;
OPLR_r_off = 0.35;

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


%average sediment dry bulk density, consisted with M1, M2
rho0 = 0.5*1000*1e6; %mg/m3

%POM proportion in the sediment
POMP = 0.005;

%proportion of OC,ON and OP in sediment, 
%assumed Redfield Ratio and 0.5% OM in the estuarine sediment.

for i=1:r
    for j=1:c
        NR = 16*14/(106*12+16*14+1*31);
        PR = 31/(106*12+16*14+1*31);
        CR = 1-NR-PR;
        OCP(i,j) = CR*POMP; %106
        ONP(i,j) = NR*POMP; %16
        OPP(i,j) = PR*POMP; %1
        SP(i,j) = 0.002; % assumed decided by diatom, C:Si ~ 3
    end
end

%proportion of G1, G2 and G3
G1P = 0.005;
G2P = 0.05;
G3P = 1-G1P-G2P;

%Sediment
CTEMP = ones(r,c)*25;

CPOPG1 = ones(r,c)*rho0.*OPP*G1P;
CPOPG2 = ones(r,c)*rho0.*OPP*G2P;
CPOPG3 = ones(r,c)*rho0.*OPP*G3P;

CPONG1 = ones(r,c)*rho0.*ONP*G1P;
CPONG2 = ones(r,c)*rho0.*ONP*G2P;
CPONG3 = ones(r,c)*rho0.*ONP*G3P;

CPOCG1 = ones(r,c)*rho0.*OCP*G1P;
CPOCG2 = ones(r,c)*rho0.*OCP*G2P;
CPOCG3 = ones(r,c)*rho0.*OCP*G3P;

CPOS = ones(r,c)*rho0.*SP;


%Layer 1 aerobic
PO4T1TM1S = ones(r,c)*18000;

NH4T1TM1S = ones(r,c)*40;

NO3T1TM1S = ones(r,c)*2;

SIT1TM1S = ones(r,c)*300000;

CH4T1TM1S = ones(r,c)*0.0;

HST1TM1S = ones(r,c)*0.01;

%Layer 2 anaerobic
PO4T2TM1S = ones(r,c)*4400;

NH4T2TM1S = ones(r,c)*400;

NO3T2TM1S = ones(r,c)*0.1;

SIT2TM1S = ones(r,c)*300000;

HST2TM1S = ones(r,c)*1400;
CH4T2TM1S = ones(r,c)*0.005;
SO4T2TM1S = ones(r,c)*1000;

BNTHSTR1S = ones(r,c)*4;

CPOP(:,:,1)  =CPOPG1;
CPOP(:,:,2)  =CPOPG2;
CPOP(:,:,3)  =CPOPG3;
CPON(:,:,1)  =CPONG1;
CPON(:,:,2)  =CPONG2;
CPON(:,:,3)  =CPONG3;
CPOC(:,:,1)  =CPOCG1;
CPOC(:,:,2)  =CPOCG2;
CPOC(:,:,3)  =CPOCG3;

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