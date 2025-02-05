clear all;
addpath(path,'../NC_file_generation/');
init_file = 'biotest_ini_bio.nc'; delete(init_file);
year = 0;
ot_start = 0.0;
grd_name =  './bio1dtest_grd.nc';
grd2 = '../Model_grid/ROMS_WFS_Piney.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N =21;
N=gn.N;
Nbed = 0;
NNS = 0;
NCS =0;
Nveg=0;
NBT=26;
bio_sed = 1;
create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,Nveg,NBT,0,bio_sed);
const_initialize(init_file,0);


[r,c] = size(lon);

Vtransform = ncread(grd2,'Vtransform');                         
Vstretching = ncread(grd2,'Vstretching');
theta_s = ncread(grd2,'theta_s');                     
theta_b = ncread(grd2,'theta_b');                     
Tcline = ncread(grd2,'Tcline');
hc = ncread(grd2,'hc');
Cs_w = ncread(grd2,'Cs_w');
Cs_r = ncread(grd2,'Cs_r');
s_w = ncread(grd2,'s_w');
s_rho = ncread(grd2,'s_rho');

spherical = ncread(grd_name,'spherical');
if(strcmp(spherical,'T'))
    ncwrite(init_file,'spherical',1);
else
    ncwrite(init_file,'spherical',0);
end

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

HSED = ones(r,c)*0.15;
VSED = ones(r,c)*0.25;
VDMIX = ones(r,c)*0.001;
VPMIX = ones(r,c)*0.00012;



ncwrite(init_file,'theta_s',theta_s);
ncwrite(init_file,'theta_b',theta_b);
ncwrite(init_file,'Tcline',Tcline);
ncwrite(init_file,'Cs_r',Cs_r);
ncwrite(init_file,'Cs_w',Cs_w);
ncwrite(init_file,'s_w',s_w);
ncwrite(init_file,'s_rho',s_rho);
ncwrite(init_file,'hc',hc);
ncwrite(init_file,'Vtransform',Vtransform);
ncwrite(init_file,'Vstretching',Vstretching);
ncwrite(init_file,'spherical',spherical);

ncwrite(init_file,'ocean_time',ot_start);
ncwrite(init_file,'salt',0*repmat(ones(size(lon)),1,1,21));
ncwrite(init_file,'temp',25*repmat(ones(size(lon)),1,1,21));
ncwrite(init_file,'SAL',0*repmat(ones(size(lon)),1,1,21));
ncwrite(init_file,'DO',7*repmat(ones(size(lon)),1,1,21));

if(bio_sed>0)
    ncwrite(init_file,'CTEMP',CTEMP);
    ncwrite(init_file,'CPOPG1',CPOPG1);
    ncwrite(init_file,'CPOPG2',CPOPG2);
    ncwrite(init_file,'CPOPG3',CPOPG3);
    ncwrite(init_file,'CPONG1',CPONG1);
    ncwrite(init_file,'CPONG2',CPONG2);
    ncwrite(init_file,'CPONG3',CPONG3);
    ncwrite(init_file,'CPOCG1',CPOCG1);
    ncwrite(init_file,'CPOCG2',CPOCG2);
    ncwrite(init_file,'CPOCG3',CPOCG3);
    ncwrite(init_file,'CPOS',CPOS);
    ncwrite(init_file,'BNTHSTR1S',BNTHSTR1S);
    ncwrite(init_file,'PO4T2TM1S',PO4T2TM1S);
    ncwrite(init_file,'NH4T2TM1S',NH4T2TM1S);
    ncwrite(init_file,'NO3T2TM1S',NO3T2TM1S);
    ncwrite(init_file,'HST2TM1S',HST2TM1S);
    ncwrite(init_file,'SIT2TM1S',SIT2TM1S);
    ncwrite(init_file,'CH4T2TM1S',CH4T2TM1S);
    ncwrite(init_file,'SO4T2TM1S',SO4T2TM1S);
    ncwrite(init_file,'PO4T1TM1S',PO4T1TM1S);
    ncwrite(init_file,'NH4T1TM1S',NH4T1TM1S);
    ncwrite(init_file,'NO3T1TM1S',NO3T1TM1S);
    ncwrite(init_file,'HST1TM1S',HST1TM1S);
    ncwrite(init_file,'SIT1TM1S',SIT1TM1S);
    ncwrite(init_file,'CH4T1TM1S',CH4T1TM1S);

    ncwrite(init_file,'HSED',HSED);
    ncwrite(init_file,'VSED',VSED);
    ncwrite(init_file,'VPMIX',VPMIX);
    ncwrite(init_file,'VDMIX',VDMIX);
end
