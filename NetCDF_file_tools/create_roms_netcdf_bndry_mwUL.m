function create_roms_netcdf_bndry_mwUL(fn,gn,t_clim,t_clim2,t_clim3,dataset_source)
%
% jcw April 18, 2009
% updated 01Sep2015 to have all 4 sides
%    for all the vars.
%

% Get some grid info. 
  [LP,MP]=size(gn.lon_rho);
  L=LP-1;
  Lm=L-1;
  M=MP-1;
  Mm=M-1;
  L  = Lm+1;
  M  = Mm+1;
  xpsi  = L;
  xrho  = LP;
  xu    = L;
  xv    = LP;
  epsi = M;
  erho = MP;
  eu   = MP;
  ev   = M;
  s    = gn.N;

%% create bndry file
nc_bndry=netcdf.create(fn,'clobber');
if isempty(nc_bndry), return, end

%% Global attributes:
myName= mfilename;
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by ',myName,' on ' datestr(now)]);
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'type', dataset_source);
%% Dimensions:

disp(' ## Defining Dimensions...')
 
psidimID = netcdf.defDim(nc_bndry,'xpsi',L);
xrhodimID = netcdf.defDim(nc_bndry,'xrho',LP);
xudimID = netcdf.defDim(nc_bndry,'xu',L);
xvdimID = netcdf.defDim(nc_bndry,'xv',LP);

epsidimID = netcdf.defDim(nc_bndry,'epsi',M);
erhodimID = netcdf.defDim(nc_bndry,'erho',MP);
eudimID = netcdf.defDim(nc_bndry,'eu',MP);
evdimID = netcdf.defDim(nc_bndry,'ev',M);
s_rhodimID = netcdf.defDim(nc_bndry,'s_rho',s);

zttdimID = netcdf.defDim(nc_bndry,'zeta_time',t_clim2);
v2tdimID = netcdf.defDim(nc_bndry,'v2d_time',t_clim2);
v3tdimID = netcdf.defDim(nc_bndry,'v3d_time',t_clim);
sltdimID = netcdf.defDim(nc_bndry,'salt_time',t_clim);
if(~isempty(t_clim3))
    dytdimID = netcdf.defDim(nc_bndry,'dye_time',t_clim3);
end
tptdimID = netcdf.defDim(nc_bndry,'temp_time',t_clim);
% mdtdimID = netcdf.defDim(nc_bndry,'mud_time',t_clim2);
 
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
ztID = netcdf.defVar(nc_bndry,'zeta_time','float',zttdimID);
netcdf.putAtt(nc_bndry,ztID,'long_name','zeta_time');
netcdf.putAtt(nc_bndry,ztID,'units','days');
netcdf.putAtt(nc_bndry,ztID,'field','zeta_time, scalar, series');

v2ID = netcdf.defVar(nc_bndry,'v2d_time','float',v2tdimID);
netcdf.putAtt(nc_bndry,v2ID,'long_name','v2d_time');
netcdf.putAtt(nc_bndry,v2ID,'units','days');
netcdf.putAtt(nc_bndry,v2ID,'field','v2d_time, scalar, series');

v3ID = netcdf.defVar(nc_bndry,'v3d_time','float',v3tdimID);
netcdf.putAtt(nc_bndry,v3ID,'long_name','v3d_time');
netcdf.putAtt(nc_bndry,v3ID,'units','days');
netcdf.putAtt(nc_bndry,v3ID,'field','v3d_time, scalar, series');

slID = netcdf.defVar(nc_bndry,'salt_time','float',sltdimID);
netcdf.putAtt(nc_bndry,slID,'long_name','salt_time');
netcdf.putAtt(nc_bndry,slID,'units','days');
netcdf.putAtt(nc_bndry,slID,'field','salt_time, scalar, series');

tpID = netcdf.defVar(nc_bndry,'temp_time','float',tptdimID);
netcdf.putAtt(nc_bndry,tpID,'long_name','temp_time');
netcdf.putAtt(nc_bndry,tpID,'units','days');
netcdf.putAtt(nc_bndry,tpID,'field','temp_time, scalar, series');

if(~isempty(t_clim3))
    dyID = netcdf.defVar(nc_bndry,'dye_time','float',dytdimID);
    netcdf.putAtt(nc_bndry,dyID,'long_name','dye_time');
    netcdf.putAtt(nc_bndry,dyID,'units','days');
    netcdf.putAtt(nc_bndry,dyID,'field','dye_time, scalar, series');
end

% mdID = netcdf.defVar(nc_bndry,'mud_time','float',mdtdimID);
% netcdf.putAtt(nc_bndry,tpID,'long_name','time for mud climatology');
% netcdf.putAtt(nc_bndry,tpID,'units','days');
% netcdf.putAtt(nc_bndry,tpID,'field','mud_time, scalar, series');

zetsID = netcdf.defVar(nc_bndry,'zeta_south','float',[xrhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zetsID,'long_name','free-surface southern boundary condition');
netcdf.putAtt(nc_bndry,zetsID,'units','meter');
netcdf.putAtt(nc_bndry,zetsID,'field','zeta_south, scalar, series');

zeteID = netcdf.defVar(nc_bndry,'zeta_east','float',[erhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zeteID,'long_name','free-surface eastern boundary condition');
netcdf.putAtt(nc_bndry,zeteID,'units','meter');
netcdf.putAtt(nc_bndry,zeteID,'field','zeta_east, scalar, series');

zetwID = netcdf.defVar(nc_bndry,'zeta_west','float',[erhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zetwID,'long_name','free-surface western boundary condition');
netcdf.putAtt(nc_bndry,zetwID,'units','meter');
netcdf.putAtt(nc_bndry,zetwID,'field','zeta_west, scalar, series');

zetnID = netcdf.defVar(nc_bndry,'zeta_north','float',[xrhodimID zttdimID]);
netcdf.putAtt(nc_bndry,zetnID,'long_name','free-surface northern boundary condition');
netcdf.putAtt(nc_bndry,zetnID,'units','meter');
netcdf.putAtt(nc_bndry,zetnID,'field','zeta_north, scalar, series');

ubsID = netcdf.defVar(nc_bndry,'ubar_south','float',[xudimID v2tdimID]);
netcdf.putAtt(nc_bndry,ubsID,'long_name','2D u-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,ubsID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ubsID,'field','ubar_south, scalar, series');

ubeID = netcdf.defVar(nc_bndry,'ubar_east','float',[eudimID v2tdimID]);
netcdf.putAtt(nc_bndry,ubeID,'long_name','2D u-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,ubeID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ubeID,'field','ubar_east, scalar, series');

ubwID = netcdf.defVar(nc_bndry,'ubar_west','float',[eudimID v2tdimID]);
netcdf.putAtt(nc_bndry,ubwID,'long_name','2D u-momentum western boundary condition');
netcdf.putAtt(nc_bndry,ubwID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ubwID,'field','ubar_west, scalar, series');

ubnID = netcdf.defVar(nc_bndry,'ubar_north','float',[xudimID v2tdimID]);
netcdf.putAtt(nc_bndry,ubnID,'long_name','2D u-momentum northern boundary condition');
netcdf.putAtt(nc_bndry,ubnID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ubnID,'field','ubar_north, scalar, series');

vbsID = netcdf.defVar(nc_bndry,'vbar_south','float',[xvdimID v2tdimID]);
netcdf.putAtt(nc_bndry,vbsID,'long_name','2D v-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,vbsID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vbsID,'field','vbar_south, scalar, series');

vbeID = netcdf.defVar(nc_bndry,'vbar_east','float',[evdimID v2tdimID]);
netcdf.putAtt(nc_bndry,vbeID,'long_name','2D v-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,vbeID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vbeID,'field','vbar_east, scalar, series');

vbwID = netcdf.defVar(nc_bndry,'vbar_west','float',[evdimID v2tdimID]);
netcdf.putAtt(nc_bndry,vbwID,'long_name','2D v-momentum western boundary condition');
netcdf.putAtt(nc_bndry,vbwID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vbwID,'field','vbar_west, scalar, series');

vbnID = netcdf.defVar(nc_bndry,'vbar_north','float',[xvdimID v2tdimID]);
netcdf.putAtt(nc_bndry,vbnID,'long_name','2D v-momentum northern boundary condition');
netcdf.putAtt(nc_bndry,vbnID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vbnID,'field','vbar_north, scalar, series');

usID = netcdf.defVar(nc_bndry,'u_south','float',[xudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,usID,'long_name','3D u-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,usID,'units','meter second-1');
netcdf.putAtt(nc_bndry,usID,'field','u_south, scalar, series');

ueID = netcdf.defVar(nc_bndry,'u_east','float',[eudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,ueID,'long_name','3D u-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,ueID,'units','meter second-1');
netcdf.putAtt(nc_bndry,ueID,'field','u_east, scalar, series');

uwID = netcdf.defVar(nc_bndry,'u_west','float',[eudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,uwID,'long_name','3D u-momentum western boundary condition');
netcdf.putAtt(nc_bndry,uwID,'units','meter second-1');
netcdf.putAtt(nc_bndry,uwID,'field','u_west, scalar, series');

unID = netcdf.defVar(nc_bndry,'u_north','float',[xudimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,unID,'long_name','3D u-momentum northern boundary condition');
netcdf.putAtt(nc_bndry,unID,'units','meter second-1');
netcdf.putAtt(nc_bndry,unID,'field','u_north, scalar, series');

vsID = netcdf.defVar(nc_bndry,'v_south','float',[xvdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,vsID,'long_name','3D v-momentum southern boundary condition');
netcdf.putAtt(nc_bndry,vsID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vsID,'field','v_south, scalar, series');

veID = netcdf.defVar(nc_bndry,'v_east','float',[evdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,veID,'long_name','3D v-momentum eastern boundary condition');
netcdf.putAtt(nc_bndry,veID,'units','meter second-1');
netcdf.putAtt(nc_bndry,veID,'field','v_east, scalar, series');

vwID = netcdf.defVar(nc_bndry,'v_west','float',[evdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,vwID,'long_name','3D v-momentum western boundary condition');
netcdf.putAtt(nc_bndry,vwID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vwID,'field','v_west, scalar, series');

vnID = netcdf.defVar(nc_bndry,'v_north','float',[xvdimID s_rhodimID v3tdimID]);
netcdf.putAtt(nc_bndry,vnID,'long_name','3D v-momentum northern boundary condition');
netcdf.putAtt(nc_bndry,vnID,'units','meter second-1');
netcdf.putAtt(nc_bndry,vnID,'field','v_north, scalar, series');

tmpsID = netcdf.defVar(nc_bndry,'temp_south','float',[xrhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc_bndry,tmpsID,'long_name','3D temperature southern boundary condition');
netcdf.putAtt(nc_bndry,tmpsID,'units','C');
netcdf.putAtt(nc_bndry,tmpsID,'field','temp_south, scalar, series');

tmpeID = netcdf.defVar(nc_bndry,'temp_east','float',[erhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc_bndry,tmpeID,'long_name','3D temperature eastern boundary condition');
netcdf.putAtt(nc_bndry,tmpeID,'units','C');
netcdf.putAtt(nc_bndry,tmpeID,'field','temp_east, scalar, series');

tmpwID = netcdf.defVar(nc_bndry,'temp_west','float',[erhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc_bndry,tmpwID,'long_name','3D temperature western boundary condition');
netcdf.putAtt(nc_bndry,tmpwID,'units','C');
netcdf.putAtt(nc_bndry,tmpwID,'field','temp_west, scalar, series');

tmpnID = netcdf.defVar(nc_bndry,'temp_north','float',[xrhodimID s_rhodimID tptdimID]);
netcdf.putAtt(nc_bndry,tmpnID,'long_name','3D temperature northern boundary condition');
netcdf.putAtt(nc_bndry,tmpnID,'units','C');
netcdf.putAtt(nc_bndry,tmpnID,'field','temp_north, scalar, series');

salsID = netcdf.defVar(nc_bndry,'salt_south','float',[xrhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc_bndry,salsID,'long_name','3D salinity southern boundary condition');
netcdf.putAtt(nc_bndry,salsID,'units','psu');
netcdf.putAtt(nc_bndry,salsID,'field','salt_south, scalar, series');

saleID = netcdf.defVar(nc_bndry,'salt_east','float',[erhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc_bndry,saleID,'long_name','3D salinity eastern boundary condition');
netcdf.putAtt(nc_bndry,saleID,'units','psu');
netcdf.putAtt(nc_bndry,saleID,'field','salt_east, scalar, series');

salwID = netcdf.defVar(nc_bndry,'salt_west','float',[erhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc_bndry,salwID,'long_name','3D salinity western boundary condition');
netcdf.putAtt(nc_bndry,salwID,'units','psu');
netcdf.putAtt(nc_bndry,salwID,'field','salt_west, scalar, series');

salnID = netcdf.defVar(nc_bndry,'salt_north','float',[xrhodimID s_rhodimID sltdimID]);
netcdf.putAtt(nc_bndry,salnID,'long_name','3D salinity northern boundary condition');
netcdf.putAtt(nc_bndry,salnID,'units','psu');
netcdf.putAtt(nc_bndry,salnID,'field','salt_north, scalar, series');

if(~isempty(t_clim3))
dye01sID = netcdf.defVar(nc_bndry,'dye_south_01','float',[xrhodimID s_rhodimID dytdimID]);
netcdf.putAtt(nc_bndry,dye01sID,'long_name','3D tracer 01 southern boundary condition');
netcdf.putAtt(nc_bndry,dye01sID,'units','mg/L');
netcdf.putAtt(nc_bndry,dye01sID,'field','dye_south_01, scalar, series');

dye01eID = netcdf.defVar(nc_bndry,'dye_east_01','float',[erhodimID s_rhodimID dytdimID]);
netcdf.putAtt(nc_bndry,dye01eID,'long_name','3D tracer 01 eastern boundary condition');
netcdf.putAtt(nc_bndry,dye01eID,'units','mg/L');
netcdf.putAtt(nc_bndry,dye01eID,'field','dye_east_01, scalar, series');

%dye01wID = netcdf.defVar(nc_bndry,'dye_west_01','float',[erhodimID s_rhodimID dytdimID]);
%netcdf.putAtt(nc_bndry,dye01wID,'long_name','3D tracer 01 western boundary condition');
%netcdf.putAtt(nc_bndry,dye01wID,'units','mg/L');
%netcdf.putAtt(nc_bndry,dye01wID,'field','dye_west_01, scalar, series');

dye01nID = netcdf.defVar(nc_bndry,'dye_north_01','float',[xrhodimID s_rhodimID dytdimID]);
netcdf.putAtt(nc_bndry,dye01nID,'long_name','3D tracer 01 northern boundary condition');
netcdf.putAtt(nc_bndry,dye01nID,'units','mg/L');
netcdf.putAtt(nc_bndry,dye01nID,'field','dye_north_01, scalar, series');
end

% mudn01ID = netcdf.defVar(nc_bndry,'mud_north_01','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudn01ID,'long_name','northern boundary mud,size calss 01');
% netcdf.putAtt(nc_bndry,mudn01ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudn01ID,'field','mud_north_01, scalar, series');
% 
% mudn02ID = netcdf.defVar(nc_bndry,'mud_north_02','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudn02ID,'long_name','northern boundary mud,size calss 02');
% netcdf.putAtt(nc_bndry,mudn02ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudn02ID,'field','mud_north_02, scalar, series');
% 
% mudn03ID = netcdf.defVar(nc_bndry,'mud_north_03','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudn03ID,'long_name','northern boundary mud,size calss 03');
% netcdf.putAtt(nc_bndry,mudn03ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudn03ID,'field','mud_north_03, scalar, series');
% 
% mudn04ID = netcdf.defVar(nc_bndry,'mud_north_04','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudn04ID,'long_name','northern boundary mud,size calss 04');
% netcdf.putAtt(nc_bndry,mudn04ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudn04ID,'field','mud_north_04, scalar, series');
% 
% mudn05ID = netcdf.defVar(nc_bndry,'mud_north_05','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudn05ID,'long_name','northern boundary mud,size calss 05');
% netcdf.putAtt(nc_bndry,mudn05ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudn05ID,'field','mud_north_05, scalar, series');
% 
% muds01ID = netcdf.defVar(nc_bndry,'mud_south_01','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,muds01ID,'long_name','southern boundary mud,size calss 01');
% netcdf.putAtt(nc_bndry,muds01ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,muds01ID,'field','mud_south_01, scalar, series');
% 
% muds02ID = netcdf.defVar(nc_bndry,'mud_south_02','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,muds02ID,'long_name','southern boundary mud,size calss 02');
% netcdf.putAtt(nc_bndry,muds02ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,muds02ID,'field','mud_south_02, scalar, series');
% 
% muds03ID = netcdf.defVar(nc_bndry,'mud_south_03','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,muds03ID,'long_name','southern boundary mud,size calss 03');
% netcdf.putAtt(nc_bndry,muds03ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,muds03ID,'field','mud_south_03, scalar, series');
% 
% muds04ID = netcdf.defVar(nc_bndry,'mud_south_04','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,muds04ID,'long_name','southern boundary mud,size calss 04');
% netcdf.putAtt(nc_bndry,muds04ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,muds04ID,'field','mud_south_04, scalar, series');
% 
% muds05ID = netcdf.defVar(nc_bndry,'mud_south_05','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,muds05ID,'long_name','southern boundary mud,size calss 05');
% netcdf.putAtt(nc_bndry,muds05ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,muds05ID,'field','mud_south_05, scalar, series');
% 
% mude01ID = netcdf.defVar(nc_bndry,'mud_east_01','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mude01ID,'long_name','eastern boundary mud,size calss 01');
% netcdf.putAtt(nc_bndry,mude01ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mude01ID,'field','mud_east_01, scalar, series');
% 
% mude02ID = netcdf.defVar(nc_bndry,'mud_east_02','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mude02ID,'long_name','eastern boundary mud,size calss 02');
% netcdf.putAtt(nc_bndry,mude02ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mude02ID,'field','mud_east_02, scalar, series');
% 
% mude03ID = netcdf.defVar(nc_bndry,'mud_east_03','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mude03ID,'long_name','eastern boundary mud,size calss 03');
% netcdf.putAtt(nc_bndry,mude03ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mude03ID,'field','mud_east_03, scalar, series');
% 
% mude04ID = netcdf.defVar(nc_bndry,'mud_east_04','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mude04ID,'long_name','eastern boundary mud,size calss 04');
% netcdf.putAtt(nc_bndry,mude04ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mude04ID,'field','mud_east_04, scalar, series');
% 
% mude05ID = netcdf.defVar(nc_bndry,'mud_east_05','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mude05ID,'long_name','eastern boundary mud,size calss 05');
% netcdf.putAtt(nc_bndry,mude05ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mude05ID,'field','mud_east_05, scalar, series');
% 
% mudw01ID = netcdf.defVar(nc_bndry,'mud_west_01','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudw01ID,'long_name','western boundary mud,size calss 01');
% netcdf.putAtt(nc_bndry,mudw01ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudw01ID,'field','mud_west_01, scalar, series');
% 
% mudw02ID = netcdf.defVar(nc_bndry,'mud_west_02','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudw02ID,'long_name','western boundary mud,size calss 02');
% netcdf.putAtt(nc_bndry,mudw02ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudw02ID,'field','mud_west_02, scalar, series');
% 
% mudw03ID = netcdf.defVar(nc_bndry,'mud_west_03','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudw03ID,'long_name','western boundary mud,size calss 03');
% netcdf.putAtt(nc_bndry,mudw03ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudw03ID,'field','mud_west_03, scalar, series');
% 
% mudw04ID = netcdf.defVar(nc_bndry,'mud_west_04','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudw04ID,'long_name','western boundary mud,size calss 04');
% netcdf.putAtt(nc_bndry,mudw04ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudw04ID,'field','mud_west_04, scalar, series');
% 
% mudw05ID = netcdf.defVar(nc_bndry,'mud_west_05','float',[xrhodimID s_rhodimID mdtdimID]);
% netcdf.putAtt(nc_bndry,mudw05ID,'long_name','western boundary mud,size calss 05');
% netcdf.putAtt(nc_bndry,mudw05ID,'units','kilogram meter-3');
% netcdf.putAtt(nc_bndry,mudw05ID,'field','mud_west_05, scalar, series');

%close file
netcdf.close(nc_bndry)

