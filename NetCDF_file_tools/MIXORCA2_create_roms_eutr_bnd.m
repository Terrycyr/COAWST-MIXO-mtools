clear all; close all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\NC_file_generation');

grd_name =  '../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
fn = 'WFS_2002_2003_bry_bio_mixo.nc';

dir_flag = [1 1 0 1]; %S N W E

date_bry = datenum(2002,10,1,0,0,0):1:datenum(2003,12,1,0,0,0);
time_ref = datenum(2002,11,1);
year = datevec(date_bry(1));
year = year(1);
bio_time = date_bry - time_ref;
t = length(date_bry);
create_roms_netcdf_bndry_eutr(fn,gn,t,dir_flag)
const_initialize(fn,0.);
nutri_adjust = 1.0;

year_end = datevec(date_bry(end));
year_end = year_end(1);

if(year_end==year+1)
    load(strcat('../GOM_preprocessing/','ocean_nutrients_bnd_',num2str(year+1),'.mat'));
    s_do1 = s_do;
    s_don1 = s_don;
    s_nh41 = s_nh4;
    s_no31 =s_no3;
    s_po41 = s_po4;
    s_sal1 = s_sal;
    s_si1 = s_si;
    s_temp1 = s_temp;

    n_do1 = n_do;
    n_don1 = n_don;
    n_nh41 = n_nh4;
    n_no31 =n_no3;
    n_po41 = n_po4;
    n_sal1 = n_sal;
    n_si1 = n_si;
    n_temp1 = n_temp;

    e_do1 = e_do;
    e_don1 = e_don;
    e_nh41 = e_nh4;
    e_no31 =e_no3;
    e_po41 = e_po4;
    e_sal1 = e_sal;
    e_si1 = e_si;
    e_temp1 = e_temp;
end

load(strcat('../GOM_preprocessing/','ocean_nutrients_bnd_',num2str(year),'.mat'));
pos_start = find(date_out==date_bry(1));
if(date_bry(end)>date_out(end))
    pos_end = length(date_out);
else
    pos_end = find(date_out==date_bry(end));
end
s_do = s_do(:,:,pos_start:pos_end);
s_don = s_don(:,:,pos_start:pos_end);
s_nh4 = s_nh4(:,:,pos_start:pos_end);
s_no3 = s_no3(:,:,pos_start:pos_end);
s_po4 = s_po4(:,:,pos_start:pos_end);
s_sal = s_sal(:,:,pos_start:pos_end);
s_si = s_si(:,:,pos_start:pos_end);
s_temp = s_temp(:,:,pos_start:pos_end);

n_do = n_do(:,:,pos_start:pos_end);
n_don = n_don(:,:,pos_start:pos_end);
n_nh4 = n_nh4(:,:,pos_start:pos_end);
n_no3 = n_no3(:,:,pos_start:pos_end);
n_po4 = n_po4(:,:,pos_start:pos_end);
n_sal = n_sal(:,:,pos_start:pos_end);
n_si = n_si(:,:,pos_start:pos_end);
n_temp = n_temp(:,:,pos_start:pos_end);

e_do = e_do(:,:,pos_start:pos_end);
e_don = e_don(:,:,pos_start:pos_end);
e_nh4 = e_nh4(:,:,pos_start:pos_end);
e_no3 = e_no3(:,:,pos_start:pos_end);
e_po4 = e_po4(:,:,pos_start:pos_end);
e_sal = e_sal(:,:,pos_start:pos_end);
e_si = e_si(:,:,pos_start:pos_end);
e_temp = e_temp(:,:,pos_start:pos_end);


if(year_end==year+1)
    merge_range1 = [size(s_do,3) size(s_do,3)+date_bry(end)-datenum(year+1,1,1)];
    merge_range2 = [size(s_do,3) size(s_do,3)+date_bry(end)-datenum(year+1,1,1)]-size(s_do,3)+1;

    s_do(:,:,merge_range1(1):merge_range1(2)) = s_do1(:,:,merge_range2(1):merge_range2(2));
    s_don(:,:,merge_range1(1):merge_range1(2)) = s_don1(:,:,merge_range2(1):merge_range2(2));
    s_nh4(:,:,merge_range1(1):merge_range1(2)) = s_nh41(:,:,merge_range2(1):merge_range2(2));
    s_no3(:,:,merge_range1(1):merge_range1(2)) = s_no31(:,:,merge_range2(1):merge_range2(2));
    s_po4(:,:,merge_range1(1):merge_range1(2)) = s_po41(:,:,merge_range2(1):merge_range2(2));
    s_sal(:,:,merge_range1(1):merge_range1(2)) = s_sal1(:,:,merge_range2(1):merge_range2(2));
    s_si(:,:,merge_range1(1):merge_range1(2)) = s_si1(:,:,merge_range2(1):merge_range2(2));
    s_temp(:,:,merge_range1(1):merge_range1(2)) = s_temp1(:,:,merge_range2(1):merge_range2(2));

    n_do(:,:,merge_range1(1):merge_range1(2)) = n_do1(:,:,merge_range2(1):merge_range2(2));
    n_don(:,:,merge_range1(1):merge_range1(2)) = n_don1(:,:,merge_range2(1):merge_range2(2));
    n_nh4(:,:,merge_range1(1):merge_range1(2)) = n_nh41(:,:,merge_range2(1):merge_range2(2));
    n_no3(:,:,merge_range1(1):merge_range1(2)) = n_no31(:,:,merge_range2(1):merge_range2(2));
    n_po4(:,:,merge_range1(1):merge_range1(2)) = n_po41(:,:,merge_range2(1):merge_range2(2));
    n_sal(:,:,merge_range1(1):merge_range1(2)) = n_sal1(:,:,merge_range2(1):merge_range2(2));
    n_si(:,:,merge_range1(1):merge_range1(2)) = n_si1(:,:,merge_range2(1):merge_range2(2));
    n_temp(:,:,merge_range1(1):merge_range1(2)) = n_temp1(:,:,merge_range2(1):merge_range2(2));

    e_do(:,:,merge_range1(1):merge_range1(2)) = e_do1(:,:,merge_range2(1):merge_range2(2));
    e_don(:,:,merge_range1(1):merge_range1(2)) = e_don1(:,:,merge_range2(1):merge_range2(2));
    e_nh4(:,:,merge_range1(1):merge_range1(2)) = e_nh41(:,:,merge_range2(1):merge_range2(2));
    e_no3(:,:,merge_range1(1):merge_range1(2)) = e_no31(:,:,merge_range2(1):merge_range2(2));
    e_po4(:,:,merge_range1(1):merge_range1(2)) = e_po41(:,:,merge_range2(1):merge_range2(2));
    e_sal(:,:,merge_range1(1):merge_range1(2)) = e_sal1(:,:,merge_range2(1):merge_range2(2));
    e_si(:,:,merge_range1(1):merge_range1(2)) = e_si1(:,:,merge_range2(1):merge_range2(2));
    e_temp(:,:,merge_range1(1):merge_range1(2)) = e_temp1(:,:,merge_range2(1):merge_range2(2));
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

%Offshore scale due to distance from estuaries
scale = 0.01;

n_on = n_don*scale/ONDP_r_off;
s_on = s_don*scale/ONDP_r_off;
e_on = e_don*scale/ONDP_r_off;

n_oc = n_on*5;
s_oc = s_on*5;
e_oc = e_on*5;

[m_Ccell,a_Ccell] = get_Ccell;

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
A_C_south = ones(size(SIT_south))*1e7*a_Ccell; %1e7 cells/l

PHYT1_north = ones(size(SIT_north))*1e-6;
PHYT2_north = ones(size(SIT_north))*1e-8;
PHYT3_north = ones(size(SIT_north))*1e-5;
A_C_north = ones(size(SIT_north))*1e7*a_Ccell; %1e7 cells/l

PHYT1_east = ones(size(SIT_east))*1e-6;
PHYT2_east = ones(size(SIT_east))*1e-8;
PHYT3_east = ones(size(SIT_east))*1e-5;
A_C_east = ones(size(SIT_east))*1e7*a_Ccell; %1e7 cells/l


ncwrite(fn,'bio_time',bio_time);
vname = {'SAL','PHYT1','PHYT2','PHYT3','RPOP','LPOP','RDOP','LDOP','PO4T'...
    ,'RPON','LPON','RDON','LDON','NH4T','NO23','BSI','SIT','RPOC','LPOC'...
    ,'RDOC','LDOC','EXDOC','TA','DIC','O2EQ','DO','A_C','A_CHLC','A_NC'...
    ,'A_PC','M_AVGCU','M_C','M_CHLC','M_FC','M_FCHLC','M_NC','M_PC','MFNC','MFPC'};
for ivar = 1:39
    eval(['ncwrite(fn,''',vname{ivar},'_time'',bio_time);']);
end
%south
ncwrite(fn,'SAL_south',SAL_south);
ncwrite(fn,'SIT_south',SIT_south*nutri_adjust);
ncwrite(fn,'DO_south',DO_south);
ncwrite(fn,'NO23_south',NO23_south*nutri_adjust);
ncwrite(fn,'NH4T_south',NH4T_south*nutri_adjust);
ncwrite(fn,'PO4T_south',PO4T_south*nutri_adjust);
ncwrite(fn,'PHYT1_south',PHYT1_south);
ncwrite(fn,'PHYT2_south',PHYT2_south);
ncwrite(fn,'PHYT3_south',PHYT3_south);
ncwrite(fn,'A_C_south',A_C_south);
ncwrite(fn,'RDON_south',RDON_south*nutri_adjust);
ncwrite(fn,'LDON_south',LDON_south*nutri_adjust);
ncwrite(fn,'RPON_south',RPON_south*nutri_adjust);
ncwrite(fn,'LPON_south',LPON_south*nutri_adjust);
ncwrite(fn,'RDOC_south',RDOC_south*nutri_adjust);
ncwrite(fn,'LDOC_south',LDOC_south*nutri_adjust);
ncwrite(fn,'RPOC_south',RPOC_south*nutri_adjust);
ncwrite(fn,'LPOC_south',LPOC_south*nutri_adjust);

% %north
% ncwrite(fn,'SAL_north',SAL_north);
% ncwrite(fn,'SIT_north',SIT_north);
% ncwrite(fn,'DO_north',DO_north);
% ncwrite(fn,'NO23_north',NO23_north);
% ncwrite(fn,'NH4T_north',NH4T_north);
% ncwrite(fn,'PO4T_north',PO4T_north);

%east
ncwrite(fn,'SAL_east',SAL_east);
ncwrite(fn,'SIT_east',SIT_east*nutri_adjust);
ncwrite(fn,'DO_east',DO_east);
ncwrite(fn,'NO23_east',NO23_east*nutri_adjust);
ncwrite(fn,'NH4T_east',NH4T_east*nutri_adjust);
ncwrite(fn,'PO4T_east',PO4T_east*nutri_adjust);
ncwrite(fn,'PHYT1_east',PHYT1_east);
ncwrite(fn,'PHYT2_east',PHYT2_east);
ncwrite(fn,'PHYT3_east',PHYT3_east);
ncwrite(fn,'A_C_east',A_C_east);








