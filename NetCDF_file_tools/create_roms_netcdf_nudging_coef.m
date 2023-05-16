function create_roms_netcdf_nudging_coef(fn,gn)

[xi eta]=size(gn.lon_rho);
%Write NetCDF file using netcdf builtins for 2010a
nc=netcdf.create(fn,'clobber');
if isempty(nc), return, end

disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by updatclim on ' datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', 'climate forcing file from http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_analysis');

% Dimensions:
disp(' ## Defining Dimensions...')
%dimid = netcdf.defDim(ncid,dimname,dimlen)
LP=xi;
MP=eta;
L=LP-1;
M=MP-1;
N=0;

psidimID = netcdf.defDim(nc,'xpsi',L);
xrhodimID = netcdf.defDim(nc,'xrho',LP);
xudimID = netcdf.defDim(nc,'xu',L);
xvdimID = netcdf.defDim(nc,'xv',LP);

epsidimID = netcdf.defDim(nc,'epsi',M);
erhodimID = netcdf.defDim(nc,'erho',MP);
eudimID = netcdf.defDim(nc,'eu',MP);
evdimID = netcdf.defDim(nc,'ev',M);
s_rhodimID = netcdf.defDim(nc,'s_rho',gn.N);


% Variables and attributes:
disp(' ## Defining Variables, and Attributes...')
%varid = netcdf.defVar(ncid,varname,xtype,dimids)
%netcdf.putAtt(ncid,varid,attrname,attrvalue)


% 
% v2ID = netcdf.defVar(nc,'v2d_time','double',v2tdimID);
% netcdf.putAtt(nc,v2ID,'long_name','v2d_time');
% netcdf.putAtt(nc,v2ID,'units','days');
% netcdf.putAtt(nc,v2ID,'field','v2d_time, scalar, series');


lonID = netcdf.defVar(nc,'lon_rho','float',[xrhodimID erhodimID]);
netcdf.putAtt(nc,lonID,'long_name','lon_rho');
netcdf.putAtt(nc,lonID,'units','degrees');
netcdf.putAtt(nc,lonID,'FillValue_',100000.);
netcdf.putAtt(nc,lonID,'missing_value',100000.);
netcdf.putAtt(nc,lonID,'field','xp, scalar, series');

latID = netcdf.defVar(nc,'lat_rho','float',[xrhodimID erhodimID]);
netcdf.putAtt(nc,latID,'long_name','lon_rho');
netcdf.putAtt(nc,latID,'units','degrees');
netcdf.putAtt(nc,latID,'FillValue_',100000.);
netcdf.putAtt(nc,latID,'missing_value',100000.);
netcdf.putAtt(nc,latID,'field','yp, scalar, series');

m2coefID = netcdf.defVar(nc,'M2_NudgeCoef','float',[xrhodimID erhodimID]);
netcdf.putAtt(nc,m2coefID,'long_name','2D momentum inverse nudging coefficients');
netcdf.putAtt(nc,m2coefID,'units','day-1');
netcdf.putAtt(nc,m2coefID,'field','M2_NudgeCoef, scalar');

m3coefID = netcdf.defVar(nc,'M3_NudgeCoef','float',[xrhodimID erhodimID s_rhodimID]);
netcdf.putAtt(nc,m3coefID,'long_name','3D momentum inverse nudging coefficients');
netcdf.putAtt(nc,m3coefID,'units','day-1');
netcdf.putAtt(nc,m3coefID,'field','M3_NudgeCoef, scalar');

tcoefID = netcdf.defVar(nc,'tracer_NudgeCoef','float',[xrhodimID erhodimID s_rhodimID]);
netcdf.putAtt(nc,tcoefID,'long_name','generic tracer inverse nudging coefficients');
netcdf.putAtt(nc,tcoefID,'units','day-1');
netcdf.putAtt(nc,tcoefID,'field','tracer_NudgeCoef, scalar');

netcdf.close(nc)

