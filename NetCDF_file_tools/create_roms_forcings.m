function create_roms_forcings(lon,lat,time,fn,varargin)
% 
% Create NetCDF file using native netcdf builtins for bulk fluxes from NAM 3-hourly output  
% 
% Usage:
% create_roms_forcings(lon,lat,time,fn,varargin)
% 
% lon: longitude array
% lat: latitiude array
% time: time array
% fn:   file name
%
% Accepts any combinations of desired parameter(s):
% (You must use these specific names on the input line)
% Uwind: surface u-wind component (m/s)
% Vwind: surface v-wind component (m/s)
% Pair: surface air pressure (mbar)
% Tair: surface air temperature (Celsius)
% Qair: surface air relative humidity (%)
% rain: rain fall rate (kg/m2/s)
% swrad: solar shortwave radiation (W/m2)
% lwrad: solar longwave radiation (W/m2)
% 
% e.g. 
% create_roms_forcings(lon_rho,lat_rho,ocean_time, 'frc_bulk.nc', 'Uwind', 'Vwind')
%                      use 'Uwind' to identify that the var Uwind in
%                      your workspace will be filled into a var Uwind
%                      in the new netcdf file.
% 
% Zafer Defne   04/30/2012
% jcwarner      06/2012 modified to use the variables, not the indices
%

nc=netcdf.create(fn,'clobber');
if isempty(nc), return, end

disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by ' mfilename ' on ' datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', 'bulk fluxes forcing file');

% Dimensions:
disp(' ## Defining Dimensions...')
%dimid = netcdf.defDim(ncid,dimname,dimlen)

[ix iy]=size(lon);
t=length(time);
lon_dimID = netcdf.defDim(nc,'xrho',ix);
lat_dimID = netcdf.defDim(nc,'yrho',iy);
t_dimID = netcdf.defDim(nc,'time',t);
% onedimID = netcdf.defDim(nc,'one',1);

% Variables and attributes:
disp(' ## Defining Variables, and Attributes...')
%varid = netcdf.defVar(ncid,varname,xtype,dimids)
%netcdf.putAtt(ncid,varid,attrname,attrvalue)

tID = netcdf.defVar(nc,'time','float',t_dimID);
netcdf.putAtt(nc,tID,'long_name','atmospheric forcing time');
netcdf.putAtt(nc,tID,'units','days');
netcdf.putAtt(nc,tID,'field','time, scalar, series');

lonID = netcdf.defVar(nc,'lon','float',[lon_dimID lat_dimID]);
netcdf.putAtt(nc,lonID,'long_name','longitude');
netcdf.putAtt(nc,lonID,'units','degrees_east');
netcdf.putAtt(nc,lonID,'field','xp, scalar, series');

latID = netcdf.defVar(nc,'lat','float',[lon_dimID lat_dimID]);
netcdf.putAtt(nc,latID,'long_name','latitude');
netcdf.putAtt(nc,latID,'units','degrees_north');
netcdf.putAtt(nc,latID,'field','yp, scalar, series');

if sum(strcmpi(varargin,'Uwind'))>0
  wt_dimID = netcdf.defDim(nc,'wind_time',t);
  wtID = netcdf.defVar(nc,'wind_time','float',wt_dimID);
  netcdf.putAtt(nc,wtID,'long_name','wind_time');
  netcdf.putAtt(nc,wtID,'units','days');
  netcdf.putAtt(nc,wtID,'field','Uwind_time, scalar, series');

  UwindID = netcdf.defVar(nc,'Uwind','float',[lon_dimID lat_dimID wt_dimID]);
  netcdf.putAtt(nc,UwindID,'long_name','surface u-wind component');
  netcdf.putAtt(nc,UwindID,'units','meter second-1');
  netcdf.putAtt(nc,UwindID,'field','Uwind, scalar, series');
  netcdf.putAtt(nc,UwindID,'coordinates','lon lat');
  netcdf.putAtt(nc,UwindID,'time','wind_time');
end

if sum(strcmpi(varargin,'Vwind'))>0
  if ~sum(strcmpi(varargin,'Uwind'))>0
    wt_dimID = netcdf.defDim(nc,'wind_time',t);
    wtID = netcdf.defVar(nc,'wind_time','float',wt_dimID);
    netcdf.putAtt(nc,wtID,'long_name','wind_time');
    netcdf.putAtt(nc,wtID,'units','days');
    netcdf.putAtt(nc,wtID,'field','Vwind_time, scalar, series');
  end

  VwindID = netcdf.defVar(nc,'Vwind','float',[lon_dimID lat_dimID wt_dimID]);
  netcdf.putAtt(nc,VwindID,'long_name','surface v-wind component');
  netcdf.putAtt(nc,VwindID,'units','meter second-1');
  netcdf.putAtt(nc,VwindID,'field','Vwind, scalar, series');
  netcdf.putAtt(nc,VwindID,'coordinates','lon lat');
  netcdf.putAtt(nc,VwindID,'time','wind_time');
end

if sum(strcmpi(varargin,'Pair'))>0
  Pat_dimID = netcdf.defDim(nc,'Pair_time',t);
  PatID = netcdf.defVar(nc,'Pair_time','float',Pat_dimID);
  netcdf.putAtt(nc,PatID,'long_name','Pair_time');
  netcdf.putAtt(nc,PatID,'units','days');
  netcdf.putAtt(nc,PatID,'field','Pair_time, scalar, series');

  PairID = netcdf.defVar(nc,'Pair','float',[lon_dimID lat_dimID Pat_dimID]);
  netcdf.putAtt(nc,PairID,'long_name','surface air pressure');
  netcdf.putAtt(nc,PairID,'units','millibar');
  netcdf.putAtt(nc,PairID,'field','Pair, scalar, series');
  netcdf.putAtt(nc,PairID,'coordinates','lon lat');
  netcdf.putAtt(nc,PairID,'time','Pair_time');
end

if sum(strcmpi(varargin,'Tair'))>0
  Tat_dimID = netcdf.defDim(nc,'Tair_time',t);
  TatID = netcdf.defVar(nc,'Tair_time','float',Tat_dimID);
  netcdf.putAtt(nc,TatID,'long_name','Tair_time');
  netcdf.putAtt(nc,TatID,'units','days');
  netcdf.putAtt(nc,TatID,'field','Tair_time, scalar, series');

  TairID = netcdf.defVar(nc,'Tair','float',[lon_dimID lat_dimID Tat_dimID]);
  netcdf.putAtt(nc,TairID,'long_name','surface air temperature');
  netcdf.putAtt(nc,TairID,'units','Celsius');
  netcdf.putAtt(nc,TairID,'field','Tair, scalar, series');
  netcdf.putAtt(nc,TairID,'coordinates','lon lat');
  netcdf.putAtt(nc,TairID,'time','Tair_time');
end

if sum(strcmpi(varargin,'Qair'))>0
  Qat_dimID = netcdf.defDim(nc,'Qair_time',t);
  QatID = netcdf.defVar(nc,'Qair_time','float',Qat_dimID);
  netcdf.putAtt(nc,QatID,'long_name','Qair_time');
  netcdf.putAtt(nc,QatID,'units','days');
  netcdf.putAtt(nc,QatID,'field','Qair_time, scalar, series');

  QairID = netcdf.defVar(nc,'Qair','float',[lon_dimID lat_dimID Qat_dimID]);
  netcdf.putAtt(nc,QairID,'long_name','surface air relative humidity');
  netcdf.putAtt(nc,QairID,'units','percentage');
  netcdf.putAtt(nc,QairID,'field','Qair, scalar, series');
  netcdf.putAtt(nc,QairID,'coordinates','lon lat');
  netcdf.putAtt(nc,QairID,'time','Qair_time');
end

if sum(strcmpi(varargin,'rain'))>0
  rt_dimID = netcdf.defDim(nc,'rain_time',t);
  rtID = netcdf.defVar(nc,'rain_time','float',rt_dimID);
  netcdf.putAtt(nc,rtID,'long_name','rain_time');
  netcdf.putAtt(nc,rtID,'units','days');
  netcdf.putAtt(nc,rtID,'field','rain_time, scalar, series');

  rainID = netcdf.defVar(nc,'rain','float',[lon_dimID lat_dimID rt_dimID]);
  netcdf.putAtt(nc,rainID,'long_name','rain fall rate');
  netcdf.putAtt(nc,rainID,'units','kilogram meter-2 second-1');
  netcdf.putAtt(nc,rainID,'field','rain, scalar, series');
  netcdf.putAtt(nc,rainID,'coordinates','lon lat');
  netcdf.putAtt(nc,rainID,'time','rain_time');
end

if sum(strcmpi(varargin,'swrad'))>0
  swrt_dimID = netcdf.defDim(nc,'swrad_time',t);
  swrtID = netcdf.defVar(nc,'swrad_time','float',swrt_dimID);
  netcdf.putAtt(nc,swrtID,'long_name','swrad_time');
  netcdf.putAtt(nc,swrtID,'units','days');
  netcdf.putAtt(nc,swrtID,'field','swrad_time, scalar, series');

  swradID = netcdf.defVar(nc,'swrad','float',[lon_dimID lat_dimID swrt_dimID]);
  netcdf.putAtt(nc,swradID,'long_name','solar shortwave radiation');
  netcdf.putAtt(nc,swradID,'units','Watts meter-2');
  netcdf.putAtt(nc,swradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,swradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,swradID,'field','swrad, scalar, series');
  netcdf.putAtt(nc,swradID,'coordinates','lon lat');
  netcdf.putAtt(nc,swradID,'time','swrad_time');
end

if sum(strcmpi(varargin,'lwrad'))>0
  lwrt_dimID = netcdf.defDim(nc,'lwrad_time',t);
  lwrtID = netcdf.defVar(nc,'lwrad_time','float',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','lwrad_time');
  netcdf.putAtt(nc,lwrtID,'units','days');
  netcdf.putAtt(nc,lwrtID,'field','lwrad_time, scalar, series');

  lwradID = netcdf.defVar(nc,'lwrad','float',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','solar longwave radiation');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'positive_value','downward flux, heating');
  netcdf.putAtt(nc,lwradID,'negative_value','upward flux, cooling');
  netcdf.putAtt(nc,lwradID,'field','lwrad, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','lwrad_time');
end

if sum(strcmpi(varargin,'lwrad_down'))>0
  lwrt_dimID = netcdf.defDim(nc,'lwrad_time',t);
  lwrtID = netcdf.defVar(nc,'lwrad_time','float',lwrt_dimID);
  netcdf.putAtt(nc,lwrtID,'long_name','lwrad_time');
  netcdf.putAtt(nc,lwrtID,'units','days');
  netcdf.putAtt(nc,lwrtID,'field','lwrad_time, scalar, series');

  lwradID = netcdf.defVar(nc,'lwrad_down','float',[lon_dimID lat_dimID lwrt_dimID]);
  netcdf.putAtt(nc,lwradID,'long_name','downwelling longwave radiation flux');
  netcdf.putAtt(nc,lwradID,'units','Watts meter-2');
  netcdf.putAtt(nc,lwradID,'field','lwrad, scalar, series');
  netcdf.putAtt(nc,lwradID,'coordinates','lon lat');
  netcdf.putAtt(nc,lwradID,'time','lwrad_time');
end

if sum(strcmpi(varargin,'dye_01_sflux'))>0
  dye01t_dimID = netcdf.defDim(nc,'dye_time',t);
  dye01tID = netcdf.defVar(nc,'dye_time','float',dye01t_dimID);
  netcdf.putAtt(nc,dye01tID,'long_name','dye_time');
  netcdf.putAtt(nc,dye01tID,'units','days');
  netcdf.putAtt(nc,dye01tID,'field','dye_time, scalar, series');

  dye01ID = netcdf.defVar(nc,'dye_01_sflux','float',[lon_dimID lat_dimID dye01t_dimID]);
  netcdf.putAtt(nc,dye01ID,'long_name','tracer 01 surface flux');
  netcdf.putAtt(nc,dye01ID,'units','mg/L');
  netcdf.putAtt(nc,dye01ID,'field','dye01, scalar, series');
  netcdf.putAtt(nc,dye01ID,'coordinates','lon lat');
  netcdf.putAtt(nc,dye01ID,'time','dye_time');
end

if sum(strcmpi(varargin,'PAR'))>0
  PARt_dimID = netcdf.defDim(nc,'PAR_time',t);
  PARtID = netcdf.defVar(nc,'PAR_time','float',PARt_dimID);
  netcdf.putAtt(nc,PARtID,'long_name','PAR_time');
  netcdf.putAtt(nc,PARtID,'units','days');
  netcdf.putAtt(nc,PARtID,'field','PAR_time, scalar, series');

  PARID = netcdf.defVar(nc,'PAR','float',[lon_dimID lat_dimID PARt_dimID]);
  netcdf.putAtt(nc,PARID,'long_name','photosynthetically active radiation');
  netcdf.putAtt(nc,PARID,'units','uE meter-2 second-1');
  netcdf.putAtt(nc,PARID,'field','PAR, scalar, series');
  netcdf.putAtt(nc,PARID,'coordinates','lon lat');
  netcdf.putAtt(nc,PARID,'time','PAR_time');
end
netcdf.close(nc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%open the file for writing
nc=netcdf.open(fn,'NC_WRITE');

%write time
ID=netcdf.inqVarID(nc,'time');
netcdf.putVar(nc,ID,time);

%write lon and lat
ID=netcdf.inqVarID(nc,'lon');
netcdf.putVar(nc,ID,lon);
ID=netcdf.inqVarID(nc,'lat');
netcdf.putVar(nc,ID,lat);

if sum(strcmpi(varargin,'Uwind'))>0
  ID=netcdf.inqVarID(nc,'wind_time');
  netcdf.putVar(nc,ID,time);
  Uwind=evalin('base','Uwind');
  ID=netcdf.inqVarID(nc,'Uwind');
  netcdf.putVar(nc,ID,Uwind);
end
if sum(strcmpi(varargin,'Vwind'))>0
  ID=netcdf.inqVarID(nc,'wind_time');
  netcdf.putVar(nc,ID,time);
  Vwind=evalin('base','Vwind');
  ID=netcdf.inqVarID(nc,'Vwind');
  netcdf.putVar(nc,ID,Vwind);
end
if sum(strcmpi(varargin,'Pair'))>0
  ID=netcdf.inqVarID(nc,'Pair_time');
  netcdf.putVar(nc,ID,time);
  Pair=evalin('base','Pair');
  ID=netcdf.inqVarID(nc,'Pair');
  netcdf.putVar(nc,ID,Pair);
end
if sum(strcmpi(varargin,'Tair'))>0
  ID=netcdf.inqVarID(nc,'Tair_time');
  netcdf.putVar(nc,ID,time);
  Tair=evalin('base','Tair');
  ID=netcdf.inqVarID(nc,'Tair');
  netcdf.putVar(nc,ID,Tair);
end
if sum(strcmpi(varargin,'Qair'))>0
  ID=netcdf.inqVarID(nc,'Qair_time');
  netcdf.putVar(nc,ID,time);
  Qair=evalin('base','Qair');
  ID=netcdf.inqVarID(nc,'Qair');
  netcdf.putVar(nc,ID,Qair);
end
if sum(strcmpi(varargin,'rain'))>0
  ID=netcdf.inqVarID(nc,'rain_time');
  netcdf.putVar(nc,ID,time);
  rain=evalin('base','rain');
  ID=netcdf.inqVarID(nc,'rain');
  netcdf.putVar(nc,ID,rain);
end
if sum(strcmpi(varargin,'swrad'))>0
  ID=netcdf.inqVarID(nc,'swrad_time');
  netcdf.putVar(nc,ID,time);
  swrad=evalin('base','swrad');
  ID=netcdf.inqVarID(nc,'swrad');
  netcdf.putVar(nc,ID,swrad);
end
if sum(strcmpi(varargin,'lwrad'))>0
  ID=netcdf.inqVarID(nc,'lwrad_time');
  netcdf.putVar(nc,ID,time);
  lwrad=evalin('base','lwrad');
  ID=netcdf.inqVarID(nc,'lwrad');
  netcdf.putVar(nc,ID,lwrad);
end
if sum(strcmpi(varargin,'lwrad_down'))>0
  ID=netcdf.inqVarID(nc,'lwrad_time');
  netcdf.putVar(nc,ID,time);
  lwrad_down=evalin('base','lwrad_down');
  ID=netcdf.inqVarID(nc,'lwrad_down');
  netcdf.putVar(nc,ID,lwrad_down);
end
if sum(strcmpi(varargin,'dye_01_sflux'))>0
  ID=netcdf.inqVarID(nc,'dye_time');
  netcdf.putVar(nc,ID,time);
  dye_01_sflux=evalin('base','dye_01_sflux');
  ID=netcdf.inqVarID(nc,'dye_01_sflux');
  netcdf.putVar(nc,ID,dye_01_sflux);
end

if sum(strcmpi(varargin,'PAR'))>0
  ID=netcdf.inqVarID(nc,'PAR_time');
  netcdf.putVar(nc,ID,time);
  PAR =evalin('base','PAR');
  ID=netcdf.inqVarID(nc,'PAR');
  netcdf.putVar(nc,ID,PAR);
end
netcdf.close(nc)


