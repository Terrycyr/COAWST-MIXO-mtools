clear all;
clearvars; %close all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\NC_file_generation');
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\m_map');
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

year = 2022;

%1 uniform; 2 offshore only; high biomass; 3 bay only;  4 no syn.; 5 no coastal 
coast = 2; 

init_file = ['WFS_',num2str(year),'_hot_bio_mixo.nc']; delete(init_file); 

dataset_source = 'Data from history file of last run';
ot_start = 0.0;
grd_name =  '../../../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
bay_mask = ncread('../../../Model_grid/ROMS_WFS_new_bay_mask.nc','mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N =length(ncread(grd_name,'Cs_r'));
N=gn.N;
Nbed = 0;
NNS = 0;
NCS =0;
Nveg=0;
NBT=39;
bio_sed = 1;
sed_flag=0;
nutri_adjust = 1.0;
shelf_adjust = 1.0;
SI_INI_scale = 1.0;
MFC_INI = 0.01;
M_CP = 40;
M_CN = 5.67;
M_CCHL = 30;
kb_bg = 0;

create_roms_netcdf_init_mw(init_file,gn,Nbed,NNS,NCS,Nveg,NBT,0,bio_sed,dataset_source);
const_initialize(init_file,0);

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


fn_his = './0901_New_v4/WFS_RCA_2022_his_00019.nc';
fn_dia = './0901_New_v4/WFS_RCA_2022_dia_00019.nc';
var_his = {'salt','temp','SAL','BSI','SIT','DO','NO23','PO4T','NH4T','PHYT1',...
    'PHYT2','ZOO1','ZOO2','LDON','RDON','LPON','RPON','LDOC','RDOC','LPOC'...
    ,'RPOC','A_C','A_NC','A_PC','A_CHLC'};

var_dia ={'CTEMP','CPOPG1','CPOPG2','CPOPG3','CPONG1','CPONG2','CPONG3'...
    ,'CPOCG1','CPOCG2','CPOCG3','CPOS','BNTHSTR1S','PO4T2TM1S','NH4T2TM1S'...
    ,'NO3T2TM1S','HST2TM1S','SIT2TM1S','CH4T2TM1S','SO4T2TM1S','PO4T1TM1S'...
    ,'NH4T1TM1S','NO3T1TM1S','HST1TM1S','SIT1TM1S','CH4T1TM1S'};

[m_Ccell,a_Ccell] = get_Ccell;

%Set initials for kb
if(coast==1)
    coast_flag = ones(r,c,N);
elseif(coast==2)
    load('selected_kb_region.mat','select_flag');
    coast_flag = repmat(select_flag,1,1,N);
    coast_flag(:,:,round(N/2):end) = 0;
elseif(coast==3)
    coast_flag = repmat(bay_mask,1,1,N);
elseif(coast==4)
    coast_flag = ones(r,c,N); 
elseif(coast==5)
    coast_flag = repmat(dep>15&bay_mask==0,1,1,N);
end

%  Assume less in areas further to Tampa Bay and Charlotte Harbour. 
scale_distance0 = get_scale_distance(lon,lat,bay_mask,mask,0.,0.01);
scale_distance0(scale_distance0<0.01) = 0;
scale_distance = scale_distance0;

M_C = ones(r,c,N);

%1 uniform; 2 offshore only; high biomass; 3 bay only;  4 no syn.; 5 no coastal 
M_C = M_C*kb_bg*m_Ccell.*repmat(scale_distance,1,1,N); % <1 cell/l

if(coast==1)
    M_C(coast_flag==1) = kb_ini; % 1e3 cells/l
elseif(coast==2)
    M_C(coast_flag==1) = 6000*m_Ccell; % 1e4 cells/l, offshore
elseif(coast==3)
    M_C(coast_flag==1) = kb_ini; % ? cells/l
elseif(coast==4)
    M_C(coast_flag==1) = kb_ini; % ? cells/l
elseif(coast==5)
    M_C(coast_flag==1) = kb_ini; % ? cells/l
end

%Sediment
%  Assume less don in areas further to Tampa Bay and Charlotte Harbour. 
scale_distance1 = get_scale_distance(lon,lat,bay_mask,mask,0.01,0.001);
%  scale_distance = min(scale_distance1,scale_distance2);
scale_distance = scale_distance1;
%  average sediment dry bulk density, consisted with M1, M2
rho0 = 0.5*1000*1e6; %mg/m3
%  POM proportion in the sediment
POMP = 0.005;
%  proportion of OC,ON and OP in sediment, 
%  assumed Redfield Ratio and 0.5% OM in the estuarine sediment.
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
G1P = 0.01;
G2P = 0.06;
G3P = 1-G1P-G2P;

load(strcat('../../../Non_tidal_component_preprocessing/HYCOM/','non_tidal_ini_',num2str(year),'.mat'));
CTEMP = temp_ini_out(:,:,1)*0.8;

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

BNTHSTR1S = ones(r,c)*12;

HSED = ones(r,c)*0.10;
VSED = ones(r,c)*0.25;
VDMIX = ones(r,c)*0.001;
VPMIX = ones(r,c)*0.00012;

for i=1:length(var_his)
    try
        tmp = ncread(fn_his,var_his{i});
    catch
        tmp = zeros(r,c,N);
    end

    for j = 1:size(tmp,3)
        tmp2 = tmp(:,:,j,end);
        tmp2(isnan(tmp2)) = 0;
        eval(strcat(var_his{i},'_his(:,:,j) = tmp2;'));
    end
end

for i=1:length(var_dia)
    tmp = ncread(fn_dia,['dia_',var_dia{i}]);
    ndim = ndims(tmp);

    if(ndim==4)
        for j = 1:size(tmp,3)
            tmp2 = tmp(:,:,j,end);
            eval(strcat('tmp3=',var_dia{i},'(:,:,j);'));
            tmp2(isnan(tmp2)) = tmp3(isnan(tmp2));
            eval(strcat(var_dia{i},'(:,:,j) = tmp2;'));
        end
    elseif(ndim==3)
        tmp2 = tmp(:,:,end);
        eval(strcat('tmp3=',var_dia{i},';'));
        tmp2(isnan(tmp2)) = tmp3(isnan(tmp2));
        eval(strcat(var_dia{i},' = tmp2;'));
    end
end

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
ncwrite(init_file,'salt',salt_his);
ncwrite(init_file,'temp',temp_his);
ncwrite(init_file,'SAL',SAL_his);
ncwrite(init_file,'BSI',BSI_his);
ncwrite(init_file,'SIT',SIT_his);
ncwrite(init_file,'DO',DO_his);
ncwrite(init_file,'NO23',NO23_his*nutri_adjust);
ncwrite(init_file,'PO4T',PO4T_his*nutri_adjust);
ncwrite(init_file,'NH4T',NH4T_his*nutri_adjust);
ncwrite(init_file,'PHYT1',PHYT1_his);
ncwrite(init_file,'PHYT2',PHYT2_his);
ncwrite(init_file,'ZOO1',ZOO1_his);
ncwrite(init_file,'ZOO2',ZOO2_his);
ncwrite(init_file,'LDON',LDON_his*nutri_adjust);
ncwrite(init_file,'RDON',RDON_his*nutri_adjust);
ncwrite(init_file,'LPON',LPON_his*nutri_adjust);
ncwrite(init_file,'RPON',RDON_his*nutri_adjust);
ncwrite(init_file,'LDOC',LDOC_his);
ncwrite(init_file,'RDOC',RDOC_his);
ncwrite(init_file,'LPOC',LPOC_his);
ncwrite(init_file,'RPOC',RPOC_his);

ncwrite(init_file,'A_C',A_C_his);
ncwrite(init_file,'A_NC',A_NC_his);
ncwrite(init_file,'A_PC',A_PC_his);
ncwrite(init_file,'A_CHLC',A_CHLC_his);

ncwrite(init_file,'M_C',M_C);
ncwrite(init_file,'M_NC',M_C/M_CN);
ncwrite(init_file,'M_PC',M_C/M_CP);
ncwrite(init_file,'M_CHLC',M_C/M_CCHL);

ncwrite(init_file,'M_FC',M_C*MFC_INI);
ncwrite(init_file,'MFNC',M_C*MFC_INI/M_CN);
ncwrite(init_file,'MFPC',M_C*MFC_INI/M_CP);
ncwrite(init_file,'M_FCHLC',M_C*MFC_INI/M_CCHL);

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
m_contourf(lon,lat,PHYT1_his(:,:,1),100,'linestyle','none');
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
m_contourf(lon,lat,M_C(:,:,1),100,'linestyle','none');
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
m_contourf(lon,lat,log10(A_C_his(:,:,1)/a_Ccell),100,'linestyle','none');
caxis([6 8]);
hold on;
m_gshhs_h('patch',[.6 .6 .6]);
hold on;
m_grid('linestyle','none','xtick',[-88.5:2:-80.5],'ytick',[24.5:2:32.5],'tickstyle','dd','fontsize',11);
colorbar;
set(gcf,'color', [1 1 1]);
xlabel('Longitude','fontsize',9);
ylabel('Latitude','fontsize',9); 