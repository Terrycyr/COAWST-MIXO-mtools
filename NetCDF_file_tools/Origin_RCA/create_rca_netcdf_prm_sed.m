function create_rca_netcdf_prm_sed(fn,gn)
% Generate a blank netcdf file containing RCA Sediment parameters.
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

%% Variables and attributes:
HSEDID = netcdf.defVar(nc_prm,'HSED','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,HSEDID,'long_name','DEPTH OF SEDIMENT LAYER');
netcdf.putAtt(nc_prm,HSEDID,'units','cm');
netcdf.putAtt(nc_prm,HSEDID,'field','HSED, scalar, series');

VDMIXID = netcdf.defVar(nc_prm,'VDMIX','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,VDMIXID,'long_name','DISSOLVED MIXING VELOCITY');
netcdf.putAtt(nc_prm,VDMIXID,'units','M^2/day');
netcdf.putAtt(nc_prm,VDMIXID,'field','VDMIX, scalar, series');

VPMIXID = netcdf.defVar(nc_prm,'VPMIX','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,VPMIXID,'long_name','PARTICULATE MIXING VELOCITY');
netcdf.putAtt(nc_prm,VPMIXID,'units','M^2/day');
netcdf.putAtt(nc_prm,VPMIXID,'field','VPMIX, scalar, series');

VSEDID = netcdf.defVar(nc_prm,'VSED','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_prm,VSEDID,'long_name','SEDIMENTATION VELOCITY');
netcdf.putAtt(nc_prm,VSEDID,'units','cm/yr');
netcdf.putAtt(nc_prm,VSEDID,'field','VSED, scalar, series');

netcdf.close(nc_prm)
end

