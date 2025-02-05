function create_rca_netcdf_prm_eutr(fn,gn,itotsft,ft,windt,ketvft,tssft,turbt)
% Generate a blank netcdf file containing RCA EUTRO parameters.
% Yuren Chen, 2021/06/21, Guangzhou
%
% Get some grid info. 
  [LP,MP]=size(gn.lon_rho); 
  x_r  = LP+2;
  y_r = MP+2; 
  s    = gn.N;
%% create bndry file
nc_prm=netcdf.create(fn,'clobber');
if isempty(nc_prm), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_prm,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by cyr on ' datestr(now)]);
netcdf.putAtt(nc_prm,netcdf.getConstant('NC_GLOBAL'),'type', 'Nutrients data from Florida Water Atlas');
%% Dimensions:

disp(' ## Defining Dimensions...')
 
xrdimID = netcdf.defDim(nc_prm,'x_r',x_r);
yrdimID = netcdf.defDim(nc_prm,'y_r',y_r);
srdimID = netcdf.defDim(nc_prm,'s_r',s);
itotsftdimID = netcdf.defDim(nc_prm,'ITOTSF_time',itotsft);
ftdimID = netcdf.defDim(nc_prm,'F_time',ft);
windtdimID = netcdf.defDim(nc_prm,'WIND_time',windt);
ketvftdimID = netcdf.defDim(nc_prm,'KETVF_time',ketvft);
tsstdimID = netcdf.defDim(nc_prm,'TSS_time',tssft);
turbtdimID = netcdf.defDim(nc_prm,'Turb_time',turbt);

%% Variables and attributes:
itotsftID = netcdf.defVar(nc_prm,'ITOTSF_time','float',itotsftdimID);
netcdf.putAtt(nc_prm,itotsftID,'long_name','ITOTSF time');
netcdf.putAtt(nc_prm,itotsftID,'units','days');
netcdf.putAtt(nc_prm,itotsftID,'field','ITOTSF TIME, scalar, series');

ftID = netcdf.defVar(nc_prm,'F_time','float',ftdimID);
netcdf.putAtt(nc_prm,ftID,'long_name','F time');
netcdf.putAtt(nc_prm,ftID,'units','days');
netcdf.putAtt(nc_prm,ftID,'field','F TIME, scalar, series');

windtID = netcdf.defVar(nc_prm,'WIND_time','float',windtdimID);
netcdf.putAtt(nc_prm,windtID,'long_name','WIND time');
netcdf.putAtt(nc_prm,windtID,'units','days');
netcdf.putAtt(nc_prm,windtID,'field','WIND TIME, scalar, series');

ketvftID = netcdf.defVar(nc_prm,'KETVF_time','float',ketvftdimID);
netcdf.putAtt(nc_prm,ketvftID,'long_name','KETVF time');
netcdf.putAtt(nc_prm,ketvftID,'units','days');
netcdf.putAtt(nc_prm,ketvftID,'field','KETVF TIME, scalar, series');

tsstID = netcdf.defVar(nc_prm,'TSS_time','float',tsstdimID);
netcdf.putAtt(nc_prm,tsstID,'long_name','TSS time');
netcdf.putAtt(nc_prm,tsstID,'units','days');
netcdf.putAtt(nc_prm,tsstID,'field','TSS TIME, scalar, series');

turbtID = netcdf.defVar(nc_prm,'Turb_time','float',turbtdimID);
netcdf.putAtt(nc_prm,turbtID,'long_name','Turb time');
netcdf.putAtt(nc_prm,turbtID,'units','days');
netcdf.putAtt(nc_prm,turbtID,'field','Turb TIME, scalar, series');

fID = netcdf.defVar(nc_prm,'F','float',[ftdimID]);
netcdf.putAtt(nc_prm,fID,'long_name','fraction of daylight');
netcdf.putAtt(nc_prm,fID,'units',' ');
netcdf.putAtt(nc_prm,fID,'field','F, scalar, series');

itotsfID = netcdf.defVar(nc_prm,'ITOTSF','float',[itotsftdimID]);
netcdf.putAtt(nc_prm,itotsfID,'long_name','total daily solar radiation');
netcdf.putAtt(nc_prm,itotsfID,'units','ly day-1');
netcdf.putAtt(nc_prm,itotsfID,'field','ITOTSF, scalar, series');

kebsID = netcdf.defVar(nc_prm,'KEBS','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,kebsID,'long_name','extinction coefficient');
netcdf.putAtt(nc_prm,kebsID,'units','meter-1');
netcdf.putAtt(nc_prm,kebsID,'field','KEBS, scalar, series');

ketvfID = netcdf.defVar(nc_prm,'KETVF','float',[ketvftdimID]);
netcdf.putAtt(nc_prm,ketvfID,'long_name','disctinction coefficient');
netcdf.putAtt(nc_prm,ketvfID,'units','meter-1');
netcdf.putAtt(nc_prm,ketvfID,'field','KETVF, scalar, series');

klID = netcdf.defVar(nc_prm,'KL','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,klID,'long_name','2D mass transfer coefficient for reaeration');
netcdf.putAtt(nc_prm,klID,'units','meter day-1');
netcdf.putAtt(nc_prm,klID,'field','KL, scalar, series');

ssldsID = netcdf.defVar(nc_prm,'SSLDS','float',[xrdimID yrdimID srdimID]);
netcdf.putAtt(nc_prm,ssldsID,'long_name','concentration of suspended solids');
netcdf.putAtt(nc_prm,ssldsID,'units','mg SS L-1');
netcdf.putAtt(nc_prm,ssldsID,'field','SSLDS, scalar, series');

tssID = netcdf.defVar(nc_prm,'TSS','float',[xrdimID yrdimID tsstdimID]);
netcdf.putAtt(nc_prm,tssID,'long_name','total suspended solids');
netcdf.putAtt(nc_prm,tssID,'units','mg SS L-1');
netcdf.putAtt(nc_prm,tssID,'field','TSS, scalar, series');

turbID = netcdf.defVar(nc_prm,'Turb','float',[xrdimID yrdimID tsstdimID]);
netcdf.putAtt(nc_prm,turbID,'long_name','Turbulence');
netcdf.putAtt(nc_prm,turbID,'units','m s-1');
netcdf.putAtt(nc_prm,turbID,'field','Turb, scalar, series');

vsnet1ID = netcdf.defVar(nc_prm,'VSNET1','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,vsnet1ID,'long_name','2D settling efficientcy from water column to the bed for algal group 1--winter diatom group');
netcdf.putAtt(nc_prm,vsnet1ID,'units','meter day-1');
netcdf.putAtt(nc_prm,vsnet1ID,'field','VSNET1, scalar, series');

vsnet2ID = netcdf.defVar(nc_prm,'VSNET2','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,vsnet2ID,'long_name','2D settling efficientcy from water column to the bed for algal group 3--Fall group');
netcdf.putAtt(nc_prm,vsnet2ID,'units','meter day-1');
netcdf.putAtt(nc_prm,vsnet2ID,'field','VSNET2, scalar, series');

vsnet3ID = netcdf.defVar(nc_prm,'VSNET3','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,vsnet3ID,'long_name','2D settling efficientcy from water column to the bed for algal group 3--Fall group');
netcdf.putAtt(nc_prm,vsnet3ID,'units','meter day-1');
netcdf.putAtt(nc_prm,vsnet3ID,'field','VSNET3, scalar, series');

vsnet4ID = netcdf.defVar(nc_prm,'VSNET4','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,vsnet4ID,'long_name','2D settling efficientcy from water column to the bed for non-living POM');
netcdf.putAtt(nc_prm,vsnet4ID,'units','meter day-1');
netcdf.putAtt(nc_prm,vsnet4ID,'field','VSNET4, scalar, series');

windID = netcdf.defVar(nc_prm,'WIND','float',[windtdimID]);
netcdf.putAtt(nc_prm,windID,'long_name','wind speed');
netcdf.putAtt(nc_prm,windID,'units','meter sec-1');
netcdf.putAtt(nc_prm,windID,'field','WIND, scalar, series');

netcdf.close(nc_prm)

