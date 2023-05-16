function create_roms_netcdf_clm_eutr(fn,gn,t_clim,var)

[xi,eta]=size(gn.lon_rho);
%Write NetCDF file using netcdf builtins for 2010a
nc=netcdf.create(fn,'clobber');
if isempty(nc), return, end

disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by updatclim on ' datestr(now)]);
netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'type', 'Nutrients data from Florida Water Atlas');

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

octdimID = netcdf.defDim(nc,'ocean_time',netcdf.getConstant('NC_UNLIMITED'));

if(strcmp(var,'NO23'))
    no23tdimID = netcdf.defDim(nc,'NO23_time',t_clim);

    % Variables and attributes:
    disp(' ## Defining Variables, and Attributes...')

    ocID = netcdf.defVar(nc,'ocean_time','double',octdimID);
    netcdf.putAtt(nc,ocID,'long_name','ocean time');
    netcdf.putAtt(nc,ocID,'units','days');
    netcdf.putAtt(nc,ocID,'field','ocean_time, scalar, series');

    ntID = netcdf.defVar(nc,'NO23_time','double',no23tdimID);
    netcdf.putAtt(nc,ntID,'long_name','NO23_time');
    netcdf.putAtt(nc,ntID,'units','days');
    netcdf.putAtt(nc,ntID,'field','NO23_time, scalar, series');

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

    no23ID = netcdf.defVar(nc,'NO23','float',[xrhodimID erhodimID s_rhodimID no23tdimID]);
    netcdf.putAtt(nc,no23ID,'long_name','Nitrate and Nitrite');
    netcdf.putAtt(nc,no23ID,'units','mg N/l');
    netcdf.putAtt(nc,no23ID,'field','Nitrate and Nitrite, scalar, series');

elseif(strcmp(var,'NH4T'))
    nh4ttdimID = netcdf.defDim(nc,'NH4T_time',t_clim);

    % Variables and attributes:
    disp(' ## Defining Variables, and Attributes...')

    ocID = netcdf.defVar(nc,'ocean_time','double',octdimID);
    netcdf.putAtt(nc,ocID,'long_name','ocean time');
    netcdf.putAtt(nc,ocID,'units','days');
    netcdf.putAtt(nc,ocID,'field','ocean_time, scalar, series');

    nhtID = netcdf.defVar(nc,'NH4T_time','double',nh4ttdimID);
    netcdf.putAtt(nc,nhtID,'long_name','NH4T_time');
    netcdf.putAtt(nc,nhtID,'units','days');
    netcdf.putAtt(nc,nhtID,'field','NH4T_time, scalar, series');

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

    nh4tID = netcdf.defVar(nc,'NH4T','float',[xrhodimID erhodimID s_rhodimID nh4ttdimID]);
    netcdf.putAtt(nc,nh4tID,'long_name','total ammonium');
    netcdf.putAtt(nc,nh4tID,'units','mg N/l');
    netcdf.putAtt(nc,nh4tID,'field','total ammonium, scalar, series');

elseif(strcmp(var,'PO4T'))
    po4ttdimID = netcdf.defDim(nc,'PO4T_time',t_clim);

    % Variables and attributes:
    disp(' ## Defining Variables, and Attributes...')

    ocID = netcdf.defVar(nc,'ocean_time','double',octdimID);
    netcdf.putAtt(nc,ocID,'long_name','ocean time');
    netcdf.putAtt(nc,ocID,'units','days');
    netcdf.putAtt(nc,ocID,'field','ocean_time, scalar, series');

    ptID = netcdf.defVar(nc,'PO4T_time','double',po4ttdimID);
    netcdf.putAtt(nc,ptID,'long_name','PO4T_time');
    netcdf.putAtt(nc,ptID,'units','days');
    netcdf.putAtt(nc,ptID,'field','PO4T_time, scalar, series');

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

    po4tID = netcdf.defVar(nc,'PO4T','float',[xrhodimID erhodimID s_rhodimID po4ttdimID]);
    netcdf.putAtt(nc,po4tID,'long_name','total phosphate');
    netcdf.putAtt(nc,po4tID,'units','mg P/l');
    netcdf.putAtt(nc,po4tID,'field','total phosphate, scalar, series');

elseif(strcmp(var,'SIT'))
    sittdimID = netcdf.defDim(nc,'SIT_time',t_clim);

    % Variables and attributes:
    disp(' ## Defining Variables, and Attributes...')

    ocID = netcdf.defVar(nc,'ocean_time','double',octdimID);
    netcdf.putAtt(nc,ocID,'long_name','ocean time');
    netcdf.putAtt(nc,ocID,'units','days');
    netcdf.putAtt(nc,ocID,'field','ocean_time, scalar, series');

    stID = netcdf.defVar(nc,'SIT_time','double',sittdimID);
    netcdf.putAtt(nc,stID,'long_name','SIT_time');
    netcdf.putAtt(nc,stID,'units','days');
    netcdf.putAtt(nc,stID,'field','SIT_time, scalar, series');

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

    sitID = netcdf.defVar(nc,'SIT','float',[xrhodimID erhodimID s_rhodimID sittdimID]);
    netcdf.putAtt(nc,sitID,'long_name','total slicate');
    netcdf.putAtt(nc,sitID,'units','mg Si/l');
    netcdf.putAtt(nc,sitID,'field','total silicate, scalar, series');

elseif(strcmp(var,'DO'))
    dotdimID = netcdf.defDim(nc,'DO_time',t_clim);

    % Variables and attributes:
    disp(' ## Defining Variables, and Attributes...')

    ocID = netcdf.defVar(nc,'ocean_time','double',octdimID);
    netcdf.putAtt(nc,ocID,'long_name','ocean time');
    netcdf.putAtt(nc,ocID,'units','days');
    netcdf.putAtt(nc,ocID,'field','ocean_time, scalar, series');

    dtID = netcdf.defVar(nc,'DO_time','double',dotdimID);
    netcdf.putAtt(nc,dtID,'long_name','DO_time');
    netcdf.putAtt(nc,dtID,'units','days');
    netcdf.putAtt(nc,dtID,'field','DO_time, scalar, series');

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

    doID = netcdf.defVar(nc,'DO','float',[xrhodimID erhodimID s_rhodimID dotdimID]);
    netcdf.putAtt(nc,doID,'long_name','dissolved oxygen');
    netcdf.putAtt(nc,doID,'units','mg O2/l');
    netcdf.putAtt(nc,doID,'field','dissolved oxygen, scalar, series');
end

    netcdf.close(nc)
end