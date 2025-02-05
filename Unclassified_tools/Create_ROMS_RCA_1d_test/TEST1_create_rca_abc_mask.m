clear all;  
grd_name =  './bio1dtest_grd.nc';
grd2 = '../Model_grid/ROMS_WFS_Piney.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
dep = ncread(grd_name,'h');
gn = struct('lon_rho',lon);
gn.N =21;

fn = 'biotest_ABCmask.nc';

nc_id=netcdf.create(fn,'clobber');
if isempty(nc_id), return, end

[r,c] = size(lon);

disp(' ## Defining Dimensions...')
 
xrdimID = netcdf.defDim(nc_id,'xrho',r);
yrdimID = netcdf.defDim(nc_id,'eho',c);


amID = netcdf.defVar(nc_id,'ABCM','float',[xrdimID yrdimID]);
netcdf.putAtt(nc_id,amID,'long_name','Ambient Bottom Current Mask');
netcdf.putAtt(nc_id,amID,'field','ABCM, scalar, series');

netcdf.close(nc_id);

zero_initialize(fn);