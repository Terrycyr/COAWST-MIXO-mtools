clear all; close all;
grd_name =  '../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
fn = 'WFS_2022_bry_bio.nc';

n_river = 11;
dir_flag = [1 1 0 1]; %S N W E

date_out = datenum(2022,1,1,0,0,0):1:datenum(2022,12,31,24,0,0);
time_ref = datenum(2022,9,1,0,0,0);
year = datevec(date_out(1));
year = year(1);
%bio_time = date_out - date_out(1);
bio_time = date_out - time_ref;
t = length(date_out);
create_roms_netcdf_bndry_eutr(fn,gn,t,dir_flag)
const_initialize(fn,0.);

load(strcat('../GOM_preprocessing/','ocean_nutrients_bnd_',num2str(year),'.mat'));

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

%Offshore scale due to distance from estuaries
scale = 0.01;

n_on = n_don*scale/ONDP_r_off;
s_on = s_don*scale/ONDP_r_off;
e_on = e_don*scale/ONDP_r_off;

n_oc = n_on*5;
s_oc = s_on*5;
e_oc = e_on*5;

%
for layer = 1:gn.N
    SAL_south(:,layer,:) = s_sal(:,end-layer+1,:);
    DO_south(:,layer,:) = s_do(:,end-layer+1,:);
    NO23_south(:,layer,:) = s_no3(:,end-layer+1,:);
    NH4T_south(:,layer,:) = s_nh4(:,end-layer+1,:);
    PO4T_south(:,layer,:) = s_po4(:,end-layer+1,:);
    SIT_south(:,layer,:) = s_si(:,end-layer+1,:);
    RPON_south(:,layer,:) = s_on(:,end-layer+1,:)*(1-ONDP_r_off)*(1-ONLR_r_off);
    LPON_south(:,layer,:) = s_on(:,end-layer+1,:)*(1-ONDP_r_off)*ONLR_r_off;
    RDON_south(:,layer,:) = s_on(:,end-layer+1,:)*ONDP_r_off*(1-ONLR_r_off);
    LDON_south(:,layer,:) = s_on(:,end-layer+1,:)*ONDP_r_off*ONLR_r_off;
    RPOC_south(:,layer,:) = s_oc(:,end-layer+1,:)*(1-OCDP_r_off)*(1-OCLR_r_off);
    LPOC_south(:,layer,:) = s_oc(:,end-layer+1,:)*(1-OCDP_r_off)*OCLR_r_off;
    RDOC_south(:,layer,:) = s_oc(:,end-layer+1,:)*OCDP_r_off*(1-OCLR_r_off);
    LDOC_south(:,layer,:) = s_oc(:,end-layer+1,:)*OCDP_r_off*OCLR_r_off;

    SAL_north(:,layer,:) = n_sal(:,end-layer+1,:);
    DO_north(:,layer,:) = n_do(:,end-layer+1,:);
    NO23_north(:,layer,:) = n_no3(:,end-layer+1,:);
    NH4T_north(:,layer,:) = n_nh4(:,end-layer+1,:);
    PO4T_north(:,layer,:) = n_po4(:,end-layer+1,:);
    SIT_north(:,layer,:) = n_si(:,end-layer+1,:);
    RPON_north(:,layer,:) = n_on(:,end-layer+1,:)*(1-ONDP_r_off)*(1-ONLR_r_off);
    LPON_north(:,layer,:) = n_on(:,end-layer+1,:)*(1-ONDP_r_off)*ONLR_r_off;
    RDON_north(:,layer,:) = n_on(:,end-layer+1,:)*ONDP_r_off*(1-ONLR_r_off);
    LDON_north(:,layer,:) = n_on(:,end-layer+1,:)*ONDP_r_off*ONLR_r_off;
    RPOC_north(:,layer,:) = n_oc(:,end-layer+1,:)*(1-OCDP_r_off)*(1-OCLR_r_off);
    LPOC_north(:,layer,:) = n_oc(:,end-layer+1,:)*(1-OCDP_r_off)*OCLR_r_off;
    RDOC_north(:,layer,:) = n_oc(:,end-layer+1,:)*OCDP_r_off*(1-OCLR_r_off);
    LDOC_north(:,layer,:) = n_oc(:,end-layer+1,:)*OCDP_r_off*OCLR_r_off;

    SAL_east(:,layer,:) = e_sal(:,end-layer+1,:);
    DO_east(:,layer,:) = e_do(:,end-layer+1,:);
    NO23_east(:,layer,:) = e_no3(:,end-layer+1,:)*0.;
    NH4T_east(:,layer,:) = e_nh4(:,end-layer+1,:)*0.;
    PO4T_east(:,layer,:) = e_po4(:,end-layer+1,:)*0.;
    SIT_east(:,layer,:) = e_si(:,end-layer+1,:)*0.;
    RPON_east(:,layer,:) = e_on(:,end-layer+1,:)*(1-ONDP_r_off)*(1-ONLR_r_off)*0;
    LPON_east(:,layer,:) = e_on(:,end-layer+1,:)*(1-ONDP_r_off)*ONLR_r_off*0;
    RDON_east(:,layer,:) = e_on(:,end-layer+1,:)*ONDP_r_off*(1-ONLR_r_off)*0;
    LDON_east(:,layer,:) = e_on(:,end-layer+1,:)*ONDP_r_off*ONLR_r_off*0;
end

% NO23_south = NO23_south/10;
% NO23_north = NO23_north/10;
% NO23_east = NO23_east/10;

PHYT1_south = ones(size(SIT_south))*1e-6;
PHYT2_south = ones(size(SIT_south))*1e-8;
PHYT3_south = ones(size(SIT_south))*1e-5;

PHYT1_north = ones(size(SIT_north))*1e-6;
PHYT2_north = ones(size(SIT_north))*1e-8;
PHYT3_north = ones(size(SIT_north))*1e-5;

PHYT1_east = ones(size(SIT_east))*1e-6;
PHYT2_east = ones(size(SIT_east))*1e-8;
PHYT3_east = ones(size(SIT_east))*1e-5;


ncwrite(fn,'bio_time',bio_time);
vname = {'SAL','PHYT1','PHYT2','PHYT3','RPOP','LPOP','RDOP','LDOP','PO4T'...
    ,'RPON','LPON','RDON','LDON','NH4T','NO23','BSI','SIT','RPOC','LPOC'...
    ,'RDOC','LDOC','EXDOC','O2EQ','DO'};
for ivar = 1:length(vname)
    eval(['ncwrite(fn,''',vname{ivar},'_time'',bio_time);']);
end
%south
ncwrite(fn,'SAL_south',SAL_south);
ncwrite(fn,'SIT_south',SIT_south);
ncwrite(fn,'DO_south',DO_south);
ncwrite(fn,'NO23_south',NO23_south);
ncwrite(fn,'NH4T_south',NH4T_south);
ncwrite(fn,'PO4T_south',PO4T_south);
ncwrite(fn,'PHYT1_south',PHYT1_south);
ncwrite(fn,'PHYT2_south',PHYT2_south);
ncwrite(fn,'PHYT3_south',PHYT3_south);
ncwrite(fn,'RDON_south',RDON_south);
ncwrite(fn,'LDON_south',LDON_south);
ncwrite(fn,'RPON_south',RPON_south);
ncwrite(fn,'LPON_south',LPON_south);
ncwrite(fn,'RDOC_south',RDOC_south);
ncwrite(fn,'LDOC_south',LDOC_south);
ncwrite(fn,'RPOC_south',RPOC_south);
ncwrite(fn,'LPOC_south',LPOC_south);

% %north
% ncwrite(fn,'SAL_north',SAL_north);
% ncwrite(fn,'SIT_north',SIT_north);
% ncwrite(fn,'DO_north',DO_north);
% ncwrite(fn,'NO23_north',NO23_north);
% ncwrite(fn,'NH4T_north',NH4T_north);
% ncwrite(fn,'PO4T_north',PO4T_north);

%east
ncwrite(fn,'SAL_east',SAL_east);
ncwrite(fn,'SIT_east',SIT_east);
ncwrite(fn,'DO_east',DO_east);
ncwrite(fn,'NO23_east',NO23_east);
ncwrite(fn,'NH4T_east',NH4T_east);
ncwrite(fn,'PO4T_east',PO4T_east);
ncwrite(fn,'PHYT1_east',PHYT1_east);
ncwrite(fn,'PHYT2_east',PHYT2_east);
ncwrite(fn,'PHYT3_east',PHYT3_east);








