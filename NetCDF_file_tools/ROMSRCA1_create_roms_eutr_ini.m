clear all;
addpath(path,'./m_map');

year=2005;
%For K. brevis initial condition
date_start = datenum(2005,5,1,0,0,0);
date_end = datenum(2005,5,5,0,0,0);

init_file = ['WFS_',num2str(year),'_ini_bio.nc']; delete(init_file); 

dataset_source = 'Data from CHNEP Water Atlas and literatures.';
ot_start = 0.0;
grd_name =  './ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
bay_mask = ncread('./ROMS_WFS_new_bay_mask.nc','mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N =length(ncread(grd_name,'Cs_r'));
N=gn.N;
Nbed = 0;
NNS = 0;
NCS =0;
Nveg=0;
NBT=24;
bio_sed = 1;
sed_flag=0;
kb_bg = 100;
create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,Nveg,NBT,0,bio_sed,dataset_source);
const_initialize(init_file,0);

kb_Ccell = 3e-7; %0.3 ng/cell, Milroy et al., 2008
syn_Ccell = 2.5e-10; % 2.5e-4 ng/cell, Heldal et al., 2003

%Set initials for kb
mat_data = './FWC_kb_dat_0006.mat';
[kb_ini0,kb_lon,kb_lat,kb_dat] = gen_kb_ini(date_start,date_end,kb_bg,grd_name,mat_data,'mean');
kb_ini0 = kb_ini0.*kb_Ccell;
kb_ini0(kb_ini0<kb_bg*kb_Ccell) = kb_bg*kb_Ccell;
kb_ini = repmat(kb_ini0,1,1,N);

[r,c] = size(lon);
N= length(ncread(grd_name,'Cs_r'));

Vtransform = ncread(grd_name,'Vtransform');                         
Vstretching = ncread(grd_name,'Vstretching');
theta_s = ncread(grd_name,'theta_s');                     
theta_b = ncread(grd_name,'theta_b');                     
Tcline = ncread(grd_name,'Tcline');
hc = ncread(grd_name,'hc');
Cs_w = ncread(grd_name,'Cs_w');
Cs_r = ncread(grd_name,'Cs_r');
s_w = ncread(grd_name,'s_w');
s_rho = ncread(grd_name,'s_rho');

%Set initials flags for K.b.
coast_flag = ones(r,c,N);


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
CRBN1 = 6.3;
CRBS1 = 3.;

CRBP2 = 32;
CRBN2 = 5.67;
CRBS2 = 1e21;

CRBP3 = 42;
CRBN3 = 5.67;
CRBS3 = 1e21;

%Assume lees don in areas further to Tampa Bay and Charlotte Harbour. 
scale_distance1 = get_scale_distance(lon,lat,bay_mask,mask,0.01,0.001);

%Assume the increase of prey availability is inverse to depth
scale_distance2 = get_scale_depth(dep,bay_mask,mask,0.01,0.8);
scale_distance = scale_distance1;

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
        OCP(i,j) = CR*POMP*scale_distance(i,j); %106
        ONP(i,j) = NR*POMP*scale_distance(i,j); %16
        OPP(i,j) = PR*POMP*scale_distance(i,j); %1
        SP(i,j) = 0.02; % assumed decided by diatom, C:Si ~ 3
    end
end

%proportion of G1, G2 and G3
G1P = 0.005;
G2P = 0.05;
G3P = 1-G1P-G2P;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
load(strcat('./bay_ini_2005_05.mat'));
load(strcat('./ocean_nutrients_ini_',num2str(year),'.mat'));
load(strcat('./tide_ini_',num2str(year),'.mat'));
load(strcat('./non_tidal_ini_',num2str(year),'.mat'));

el_final = el_ini_out+tide_el_ini;
u2d_final = u2d_ini_out+tide_u_ini;
v2d_final = v2d_ini_out+tide_v_ini;
u_final = u_ini_out;
v_final = v_ini_out;

el_final(mask==0) = 0;
for i=1:N
    tmp = temp_ini_out(:,:,i);
    tmp(mask==0) =0;
    temp_ini_out(:,:,i) = tmp;
    
    tmp = s_ini_out(:,:,i);
    tmp(mask==0) =0;
    s_ini_out(:,:,i) = tmp;  
end

%
SAL_roms = s_ini_out;
DO = do;
NO23 = no3;
NH4T = nh4;
PO4T = po4;
SIT = si;
ON = don*1./ONDP_r;
OC = ON*5;
LDON = ON*ONDP_r_off*ONLR_r_off;
RDON = ON*ONDP_r_off*(1-ONLR_r_off);
LPON = ON*(1-ONDP_r_off)*ONLR_r_off;
RPON = ON*(1-ONDP_r_off)*(1-ONLR_r_off);
LDOC = OC*OCDP_r_off*OCLR_r_off;
RDOC = OC*OCDP_r_off*(1-OCLR_r_off);
LPOC = OC*(1-OCDP_r_off)*OCLR_r_off;
RPOC = OC*(1-OCDP_r_off)*(1-OCLR_r_off);

flag = 1;
PHYT1_roms = ones(size(si))*1e-4.*flag;
PHYT2_roms = kb_ini;
PHYT3_roms = ones(size(si))*1e5*syn_Ccell.*flag; 

k = 0;
for i=gn.N:-1:1
    k = k+1;

    DO_roms(:,:,k) = DO(:,:,i);

    tmp = SAL_roms(:,:,k);
    tmp2 = sal_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    SAL_roms(:,:,k) = tmp;

    tmp = NO23(:,:,i);
    tmp2 = no23_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    NO23_roms(:,:,k) = tmp;
    
    tmp = PO4T(:,:,i);
    tmp2 = po4_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    PO4T_roms(:,:,k) = tmp;
    
    tmp = SIT(:,:,i);
    tmp2 = sit_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    SIT_roms(:,:,k) = tmp;
    
    tmp = NH4T(:,:,i);
    tmp2 = nh4_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2));
    NH4T_roms(:,:,k) = tmp;

    tmp = LDON(:,:,i);
    tmp2 = don_bay(:,:,i);

    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*ONLR_r;
    LDON_roms(:,:,k) = tmp;

    tmp = RDON(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONLR_r);
    RDON_roms(:,:,k) = tmp;

    tmp = LPON(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)/ONDP_r*ONLR_r;
    LPON_roms(:,:,k) = tmp;

    tmp = RPON(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)/ONDP_r*(1-ONLR_r);
    RPON_roms(:,:,k) = tmp;

    tmp = LDOC(:,:,i);
    tmp2 = oc_bay(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*ONDP_r*ONLR_r;
    LDOC_roms(:,:,k) = tmp;

    tmp = RDOC(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*ONDP_r*(1-ONLR_r);
    RDOC_roms(:,:,k) = tmp;

    tmp = LPOC(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)*ONLR_r;
    LPOC_roms(:,:,k) = tmp;

    tmp = RPOC(:,:,i);
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*(1-ONDP_r)*(1-ONLR_r);
    RPOC_roms(:,:,k) = tmp;

    tmp = PHYT1_roms(:,:,k);
    tmp2 = chla_bay(:,:,i)*50;
    tmp(~isnan(tmp2)) = tmp2(~isnan(tmp2))*0.8;
    PHYT1_roms(:,:,k) = tmp;
end

BSI_roms = PHYT1_roms*(1/CRBS1);
SIT_roms = SIT_roms+PHYT1_roms*(1/CRBS1)+PHYT2_roms*(1/CRBS2)+PHYT3_roms*(1/CRBS3);
PO4T_roms = PO4T_roms+PHYT1_roms*(1/CRBP1)+PHYT2_roms*(1/CRBP2)+PHYT3_roms*(1/CRBP3);
NH4T_roms = NH4T_roms+PHYT1_roms*(1/CRBN1)+PHYT2_roms*(1/CRBN2)+PHYT3_roms*(1/CRBN3);

%Sediment
CTEMP = temp_ini_out(:,:,1)*0.9;

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
PO4T1TM1S(bay_mask==0) = 3600;

NH4T1TM1S = ones(r,c)*40;

NO3T1TM1S = ones(r,c)*2;

SIT1TM1S = ones(r,c)*300000;
SIT1TM1S(bay_mask==0) = 30000;

CH4T1TM1S = ones(r,c)*0.0;
HST1TM1S = ones(r,c)*0.01;

%Layer 2 anaerobic
PO4T2TM1S = ones(r,c)*4400;
PO4T2TM1S(bay_mask==0) = 1100;

NH4T2TM1S = ones(r,c)*400;

NO3T2TM1S = ones(r,c)*0.1;

SIT2TM1S = ones(r,c)*300000;
SIT2TM1S(bay_mask==0) = 30000;

HST2TM1S = ones(r,c)*1400;
CH4T2TM1S = ones(r,c)*0.005;
SO4T2TM1S = ones(r,c)*1000;

BNTHSTR1S = ones(r,c)*4.*scale_distance;

HSED = ones(r,c)*0.1;
VSED = ones(r,c)*0.25;
VDMIX = ones(r,c)*0.001;
VPMIX = ones(r,c)*0.00012;

spherical = ncread(grd_name,'spherical');
if(strcmp(spherical,'T'))
    ncwrite(init_file,'spherical',1);
else
    ncwrite(init_file,'spherical',0);
end

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
ncwrite(init_file,'salt',s_ini_out);
ncwrite(init_file,'temp',temp_ini_out);
ncwrite(init_file,'u',u_final);
ncwrite(init_file,'ubar',u2d_final);
ncwrite(init_file,'v',v_final);
ncwrite(init_file,'vbar',v2d_final);
ncwrite(init_file,'zeta',el_final);
ncwrite(init_file,'SAL',SAL_roms);
ncwrite(init_file,'BSI',BSI_roms);
ncwrite(init_file,'SIT',SIT_roms);
ncwrite(init_file,'DO',DO_roms);
ncwrite(init_file,'NO23',NO23_roms);
ncwrite(init_file,'PO4T',PO4T_roms);
ncwrite(init_file,'NH4T',NH4T_roms);
ncwrite(init_file,'PHYT1',PHYT1_roms);
ncwrite(init_file,'PHYT2',PHYT2_roms);
ncwrite(init_file,'PHYT3',PHYT3_roms);
ncwrite(init_file,'LDON',LDON_roms);
ncwrite(init_file,'RDON',RDON_roms);
ncwrite(init_file,'LPON',LPON_roms);
ncwrite(init_file,'RPON',RDON_roms);
ncwrite(init_file,'LDOC',LDOC_roms);
ncwrite(init_file,'RDOC',RDOC_roms);
ncwrite(init_file,'LPOC',LPOC_roms);
ncwrite(init_file,'RPOC',RPOC_roms);

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

figure('units','pixels','position',[500 150 1200 410]);
subplot(1,3,1);
m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
m_contourf(lon,lat,PHYT1_roms(:,:,end),100,'linestyle','none');
hold on;
m_gshhs_h('patch',[.6 .6 .6]);
hold on;
m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
colorbar;
set(gcf,'color', [1 1 1]);
xlabel('Longitude','fontsize',9);
ylabel('Latitude','fontsize',9); 

subplot(1,3,2);
m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
m_contourf(lon,lat,PHYT2_roms(:,:,end),100,'linestyle','none');
hold on;
m_gshhs_h('patch',[.6 .6 .6]);
hold on;
m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
colorbar;
set(gcf,'color', [1 1 1]);
xlabel('Longitude','fontsize',9);
ylabel('Latitude','fontsize',9); 

subplot(1,3,3);
m_proj('Mercator','lat',[24 32],'long',[-87.8 -80.2]);
m_contourf(lon,lat,PHYT3_roms(:,:,end),100,'linestyle','none');
hold on;
m_gshhs_h('patch',[.6 .6 .6]);
hold on;
m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
colorbar;
set(gcf,'color', [1 1 1]);
xlabel('Longitude','fontsize',9);
ylabel('Latitude','fontsize',9); 