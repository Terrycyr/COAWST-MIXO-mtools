function create_rca_netcdf_bndry_eutr(fn,gn,t,river,dir_flag)
% Generate a blank netcdf file containing RCA EUTRO boundary conditions.
% Yuren Chen, 2021/06/21, Guangzhou
%

% Get some grid info. 
  [LP,MP]=size(gn.lon_rho); 
  x_r  = LP+2;
  y_r = MP+2; 
  %x_r = LP;
  %y_r = MP;
  s    = gn.N;
  
%% create bndry file
nc_bndry=netcdf.create(fn,'clobber');
if isempty(nc_bndry), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by cyr on ' datestr(now)]);
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'type', 'Nutrients data from Florida Water Atlas');
%% Dimensions:

disp(' ## Defining Dimensions...')
 
xrdimID = netcdf.defDim(nc_bndry,'x_r',x_r);
yrdimID = netcdf.defDim(nc_bndry,'y_r',y_r);
srdimID = netcdf.defDim(nc_bndry,'s_r',s);
tdimID = netcdf.defDim(nc_bndry,'time',t);
opdimID = netcdf.defDim(nc_bndry,'option',1);
rdimID = netcdf.defDim(nc_bndry,'river',river);
%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

IBCOPTID = netcdf.defVar(nc_bndry,'IBCOPT','NC_INT',opdimID);
netcdf.putAtt(nc_bndry,IBCOPTID,'long_name','BC option, 1=sigma level-constant, 2=sigma level-time-varying, 3=standard level-constant, 4=standard level-time-varying');

IBCPWLOPTID = netcdf.defVar(nc_bndry,'IBCPWLOPT','NC_INT',opdimID);
netcdf.putAtt(nc_bndry,IBCPWLOPTID,'long_name','piecewise linear boundary concentration option, 0=step function approximation, 1=piecewise linear approximation');
%
tID = netcdf.defVar(nc_bndry,'bry_time','float',tdimID);
netcdf.putAtt(nc_bndry,tID,'long_name','bry_time');
netcdf.putAtt(nc_bndry,tID,'units','days');
netcdf.putAtt(nc_bndry,tID,'field','bry_time, scalar, series');

%South
if(dir_flag(1)==1)
    acsouthID = netcdf.defVar(nc_bndry,'A_C_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,acsouthID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,acsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,acsouthID,'field','A_C_south, scalar, series');

    achlcsouthID = netcdf.defVar(nc_bndry,'A_CHLC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,achlcsouthID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlcsouthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlcsouthID,'field','A_CHLC_south, scalar, series');

    ancsouthID = netcdf.defVar(nc_bndry,'A_NC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ancsouthID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,ancsouthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,ancsouthID,'field','A_NC_south, scalar, series');

    apcsouthID = netcdf.defVar(nc_bndry,'A_PC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,apcsouthID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apcsouthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apcsouthID,'field','A_PC_south, scalar, series');

    bsisouthID = netcdf.defVar(nc_bndry,'BSI_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,bsisouthID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsisouthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsisouthID,'field','bsi_south, scalar, series');

    dosouthID = netcdf.defVar(nc_bndry,'DO_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,dosouthID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,dosouthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,dosouthID,'field','do_south, scalar, series');

    exdocsouthID = netcdf.defVar(nc_bndry,'EXDOC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,exdocsouthID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdocsouthID,'field','exdoc_south, scalar, series');

    ldocsouthID = netcdf.defVar(nc_bndry,'LDOC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldocsouthID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldocsouthID,'field','ldoc_south, scalar, series');

    ldonsouthID = netcdf.defVar(nc_bndry,'LDON_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldonsouthID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldonsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldonsouthID,'field','ldon_south, scalar, series');

    ldopsouthID = netcdf.defVar(nc_bndry,'LDOP_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldopsouthID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopsouthID,'field','ldop_south, scalar, series');


    lpocsouthID = netcdf.defVar(nc_bndry,'LPOC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpocsouthID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpocsouthID,'field','lpoc_south, scalar, series');

    lponsouthID = netcdf.defVar(nc_bndry,'LPON_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lponsouthID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponsouthID,'field','lpon_south, scalar, series');

    lpopsouthID = netcdf.defVar(nc_bndry,'LPOP_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpopsouthID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopsouthID,'field','lpop_south, scalar, series');

    mavgusouthID = netcdf.defVar(nc_bndry,'M_AVGU_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mavgusouthID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavgusouthID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavgusouthID,'field','mavgu_south, scalar, series');

    mcsouthID = netcdf.defVar(nc_bndry,'M_C_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mcsouthID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mcsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mcsouthID,'field','mc_south, scalar, series');

    mchlcsouthID = netcdf.defVar(nc_bndry,'M_CHLC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mchlcsouthID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlcsouthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlcsouthID,'field','mchlc_south, scalar, series');

    mfcsouthID = netcdf.defVar(nc_bndry,'M_FC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfcsouthID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfcsouthID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfcsouthID,'field','mfc_south, scalar, series');

    mfchlcsouthID = netcdf.defVar(nc_bndry,'M_FCHLC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfchlcsouthID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlcsouthID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlcsouthID,'field','mfchlc_south, scalar, series');

    mncsouthID = netcdf.defVar(nc_bndry,'M_NC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mncsouthID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mncsouthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mncsouthID,'field','mnc_south, scalar, series');

    mpcsouthID = netcdf.defVar(nc_bndry,'M_PC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mpcsouthID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpcsouthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpcsouthID,'field','mpc_south, scalar, series');

    mfncsouthID = netcdf.defVar(nc_bndry,'M_FNC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfncsouthID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfncsouthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfncsouthID,'field','mfnc_south, scalar, series');

    mfpcsouthID = netcdf.defVar(nc_bndry,'M_FPC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfpcsouthID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpcsouthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpcsouthID,'field','mfpc_south, scalar, series');

    nh4tsouthID = netcdf.defVar(nc_bndry,'NH4T_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,nh4tsouthID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4tsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4tsouthID,'field','nh4t_south, scalar, series');

    no23southID = netcdf.defVar(nc_bndry,'NO23_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,no23southID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23southID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23southID,'field','no23_south, scalar, series')

    o2eqsouthID = netcdf.defVar(nc_bndry,'O2EQ_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,o2eqsouthID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqsouthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqsouthID,'field','o2eq_south, scalar, series')

    phyt1southID = netcdf.defVar(nc_bndry,'PHYT1_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt1southID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1southID,'field','phyt1_south, scalar, series')

    phyt2southID = netcdf.defVar(nc_bndry,'PHYT2_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt2southID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2southID,'field','phyt2_south, scalar, series')

    phyt3southID = netcdf.defVar(nc_bndry,'PHYT3_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt3southID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3southID,'field','phyt3_south, scalar, series')

    po4tsouthID = netcdf.defVar(nc_bndry,'PO4T_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,po4tsouthID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4tsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4tsouthID,'field','po4t_south, scalar, series')

    rdocsouthID = netcdf.defVar(nc_bndry,'RDOC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdocsouthID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdocsouthID,'field','rdoc_south, scalar, series')

    rdonsouthID = netcdf.defVar(nc_bndry,'RDON_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdonsouthID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdonsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdonsouthID,'field','rdon_south, scalar, series')

    rdopsouthID = netcdf.defVar(nc_bndry,'RDOP_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdopsouthID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rdopsouthID,'field','rdop_south, scalar, series')

    redocsouthID = netcdf.defVar(nc_bndry,'REDOC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,redocsouthID,'long_name','REACTIVE DISSOLVED ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,redocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,redocsouthID,'field','redoc_south, scalar, series')

    repocsouthID = netcdf.defVar(nc_bndry,'REPOC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,repocsouthID,'long_name','REACTIVE PARTICULATE ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,repocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,repocsouthID,'field','repoc_south, scalar, series')
    
    rpocsouthID = netcdf.defVar(nc_bndry,'RPOC_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpocsouthID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpocsouthID,'field','rpoc_south, scalar, series')
    
    rponsouthID = netcdf.defVar(nc_bndry,'RPON_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rponsouthID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponsouthID,'field','rpon_south, scalar, series')
    
    rpopsouthID = netcdf.defVar(nc_bndry,'RPOP_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpopsouthID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopsouthID,'field','rpop_south, scalar, series')
    
    salsouthID = netcdf.defVar(nc_bndry,'SAL_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,salsouthID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,salsouthID,'units','psu');
    netcdf.putAtt(nc_bndry,salsouthID,'field','sal_south, scalar, series')
    
    sitsouthID = netcdf.defVar(nc_bndry,'SIT_south','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,sitsouthID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,sitsouthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,sitsouthID,'field','sit_south, scalar, series')
end

%NORTH
if(dir_flag(2)==1)
    acnorthID = netcdf.defVar(nc_bndry,'A_C_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,acnorthID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,acnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,acnorthID,'field','A_C_north, scalar, series');

    achlcnorthID = netcdf.defVar(nc_bndry,'A_CHLC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,achlcnorthID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlcnorthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlcnorthID,'field','A_CHLC_north, scalar, series');

    ancnorthID = netcdf.defVar(nc_bndry,'A_NC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ancnorthID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,ancnorthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,ancnorthID,'field','A_NC_north, scalar, series');

    apcnorthID = netcdf.defVar(nc_bndry,'A_PC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,apcnorthID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apcnorthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apcnorthID,'field','A_PC_north, scalar, series');

    bsinorthID = netcdf.defVar(nc_bndry,'BSI_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,bsinorthID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsinorthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsinorthID,'field','bsi_north, scalar, series');

    donorthID = netcdf.defVar(nc_bndry,'DO_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,donorthID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,donorthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,donorthID,'field','do_north, scalar, series');

    exdocnorthID = netcdf.defVar(nc_bndry,'EXDOC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,exdocnorthID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdocnorthID,'field','exdoc_north, scalar, series');

    ldocnorthID = netcdf.defVar(nc_bndry,'LDOC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldocnorthID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldocnorthID,'field','ldoc_north, scalar, series');

    ldonnorthID = netcdf.defVar(nc_bndry,'LDON_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldonnorthID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldonnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldonnorthID,'field','ldon_north, scalar, series');

    ldopnorthID = netcdf.defVar(nc_bndry,'LDOP_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldopnorthID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopnorthID,'field','ldop_north, scalar, series');

    lpocnorthID = netcdf.defVar(nc_bndry,'LPOC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpocnorthID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpocnorthID,'field','lpoc_north, scalar, series');

    lponnorthID = netcdf.defVar(nc_bndry,'LPON_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lponnorthID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponnorthID,'field','lpon_north, scalar, series');

    lpopnorthID = netcdf.defVar(nc_bndry,'LPOP_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpopnorthID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopnorthID,'field','lpop_north, scalar, series');

    mavgunorthID = netcdf.defVar(nc_bndry,'M_AVGU_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mavgunorthID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavgunorthID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavgunorthID,'field','mavgu_north, scalar, series');

    mcnorthID = netcdf.defVar(nc_bndry,'M_C_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mcnorthID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mcnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mcnorthID,'field','mc_north, scalar, series');

    mchlcnorthID = netcdf.defVar(nc_bndry,'M_CHLC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mchlcnorthID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlcnorthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlcnorthID,'field','mchlc_north, scalar, series');

    mfcnorthID = netcdf.defVar(nc_bndry,'M_FC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfcnorthID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfcnorthID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfcnorthID,'field','mfc_north, scalar, series');

    mfchlcnorthID = netcdf.defVar(nc_bndry,'M_FCHLC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfchlcnorthID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlcnorthID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlcnorthID,'field','mfchlc_north, scalar, series');

    mncnorthID = netcdf.defVar(nc_bndry,'M_NC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mncnorthID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mncnorthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mncnorthID,'field','mnc_north, scalar, series');

    mpcnorthID = netcdf.defVar(nc_bndry,'M_PC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mpcnorthID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpcnorthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpcnorthID,'field','mpc_north, scalar, series');

    mfncnorthID = netcdf.defVar(nc_bndry,'M_FNC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfncnorthID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfncnorthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfncnorthID,'field','mfnc_north, scalar, series');

    mfpcnorthID = netcdf.defVar(nc_bndry,'M_FPC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfpcnorthID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpcnorthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpcnorthID,'field','mfpc_north, scalar, series');

    nh4tnorthID = netcdf.defVar(nc_bndry,'NH4T_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,nh4tnorthID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4tnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4tnorthID,'field','nh4t_north, scalar, series');

    no23northID = netcdf.defVar(nc_bndry,'NO23_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,no23northID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23northID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23northID,'field','no23_north, scalar, series')

    o2eqnorthID = netcdf.defVar(nc_bndry,'O2EQ_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,o2eqnorthID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqnorthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqnorthID,'field','o2eq_north, scalar, series')

    phyt1northID = netcdf.defVar(nc_bndry,'PHYT1_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt1northID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1northID,'field','phyt1_north, scalar, series')

    phyt2northID = netcdf.defVar(nc_bndry,'PHYT2_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt2northID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2northID,'field','phyt2_north, scalar, series')

    phyt3northID = netcdf.defVar(nc_bndry,'PHYT3_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt3northID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3northID,'field','phyt3_north, scalar, series')

    po4tnorthID = netcdf.defVar(nc_bndry,'PO4T_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,po4tnorthID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4tnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4tnorthID,'field','po4t_north, scalar, series')

    rdocnorthID = netcdf.defVar(nc_bndry,'RDOC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdocnorthID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdocnorthID,'field','rdoc_north, scalar, series')

    rdonnorthID = netcdf.defVar(nc_bndry,'RDON_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdonnorthID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdonnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdonnorthID,'field','rdon_north, scalar, series')

    rdopnorthID = netcdf.defVar(nc_bndry,'RDOP_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdopnorthID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rdopnorthID,'field','rdop_north, scalar, series')

    redocnorthID = netcdf.defVar(nc_bndry,'REDOC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,redocnorthID,'long_name','REACTIVE DISSOLVED ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,redocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,redocnorthID,'field','redoc_north, scalar, series')

    repocnorthID = netcdf.defVar(nc_bndry,'REPOC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,repocnorthID,'long_name','REACTIVE PARTICULATE ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,repocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,repocnorthID,'field','repoc_north, scalar, series')
    
    rpocnorthID = netcdf.defVar(nc_bndry,'RPOC_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpocnorthID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpocnorthID,'field','rpoc_north, scalar, series')
    
    rponnorthID = netcdf.defVar(nc_bndry,'RPON_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rponnorthID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponnorthID,'field','rpon_north, scalar, series')
    
    rpopnorthID = netcdf.defVar(nc_bndry,'RPOP_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpopnorthID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopnorthID,'field','rpop_north, scalar, series')
    
    salnorthID = netcdf.defVar(nc_bndry,'SAL_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,salnorthID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,salnorthID,'units','psu');
    netcdf.putAtt(nc_bndry,salnorthID,'field','sal_north, scalar, series')
    
    sitnorthID = netcdf.defVar(nc_bndry,'SIT_north','float',[xrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,sitnorthID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,sitnorthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,sitnorthID,'field','sit_north, scalar, series')
end
%WEST
if(dir_flag(3)==1)
    acwestID = netcdf.defVar(nc_bndry,'A_C_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,acwestID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,acwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,acwestID,'field','A_C_west, scalar, series');

    achlcwestID = netcdf.defVar(nc_bndry,'A_CHLC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,achlcwestID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlcwestID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlcwestID,'field','A_CHLC_west, scalar, series');

    ancwestID = netcdf.defVar(nc_bndry,'A_NC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ancwestID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,ancwestID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,ancwestID,'field','A_NC_west, scalar, series');

    apcwestID = netcdf.defVar(nc_bndry,'A_PC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,apcwestID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apcwestID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apcwestID,'field','A_PC_west, scalar, series');

    bsiwestID = netcdf.defVar(nc_bndry,'BSI_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,bsiwestID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsiwestID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsiwestID,'field','bsi_west, scalar, series');

    dowestID = netcdf.defVar(nc_bndry,'DO_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,dowestID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,dowestID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,dowestID,'field','do_west, scalar, series');

    exdocwestID = netcdf.defVar(nc_bndry,'EXDOC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,exdocwestID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdocwestID,'field','exdoc_west, scalar, series');

    ldocwestID = netcdf.defVar(nc_bndry,'LDOC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldocwestID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldocwestID,'field','ldoc_west, scalar, series');

    ldonwestID = netcdf.defVar(nc_bndry,'LDON_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldonwestID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldonwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldonwestID,'field','ldon_west, scalar, series');

    ldopwestID = netcdf.defVar(nc_bndry,'LDOP_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldopwestID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopwestID,'field','ldop_west, scalar, series');

    lpocwestID = netcdf.defVar(nc_bndry,'LPOC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpocwestID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpocwestID,'field','lpoc_west, scalar, series');

    lponwestID = netcdf.defVar(nc_bndry,'LPON_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lponwestID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponwestID,'field','lpon_west, scalar, series');

    lpopwestID = netcdf.defVar(nc_bndry,'LPOP_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpopwestID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopwestID,'field','lpop_west, scalar, series');

    mavguwestID = netcdf.defVar(nc_bndry,'M_AVGU_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mavguwestID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavguwestID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavguwestID,'field','mavgu_west, scalar, series');

    mcwestID = netcdf.defVar(nc_bndry,'M_C_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mcwestID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mcwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mcwestID,'field','mc_west, scalar, series');

    mchlcwestID = netcdf.defVar(nc_bndry,'M_CHLC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mchlcwestID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlcwestID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlcwestID,'field','mchlc_west, scalar, series');

    mfcwestID = netcdf.defVar(nc_bndry,'M_FC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfcwestID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfcwestID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfcwestID,'field','mfc_west, scalar, series');

    mfchlcwestID = netcdf.defVar(nc_bndry,'M_FCHLC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfchlcwestID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlcwestID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlcwestID,'field','mfchlc_west, scalar, series');

    mncwestID = netcdf.defVar(nc_bndry,'M_NC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mncwestID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mncwestID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mncwestID,'field','mnc_west, scalar, series');

    mpcwestID = netcdf.defVar(nc_bndry,'M_PC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mpcwestID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpcwestID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpcwestID,'field','mpc_west, scalar, series');

    mfncwestID = netcdf.defVar(nc_bndry,'M_FNC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfncwestID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfncwestID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfncwestID,'field','mfnc_west, scalar, series');

    mfpcwestID = netcdf.defVar(nc_bndry,'M_FPC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfpcwestID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpcwestID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpcwestID,'field','mfpc_west, scalar, series');

    nh4twestID = netcdf.defVar(nc_bndry,'NH4T_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,nh4twestID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4twestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4twestID,'field','nh4t_west, scalar, series');

    no23westID = netcdf.defVar(nc_bndry,'NO23_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,no23westID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23westID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23westID,'field','no23_west, scalar, series')

    o2eqwestID = netcdf.defVar(nc_bndry,'O2EQ_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,o2eqwestID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqwestID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqwestID,'field','o2eq_west, scalar, series')

    phyt1westID = netcdf.defVar(nc_bndry,'PHYT1_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt1westID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1westID,'field','phyt1_west, scalar, series')

    phyt2westID = netcdf.defVar(nc_bndry,'PHYT2_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt2westID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2westID,'field','phyt2_west, scalar, series')

    phyt3westID = netcdf.defVar(nc_bndry,'PHYT3_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt3westID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3westID,'field','phyt3_west, scalar, series')

    po4twestID = netcdf.defVar(nc_bndry,'PO4T_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,po4twestID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4twestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4twestID,'field','po4t_west, scalar, series')

    rdocwestID = netcdf.defVar(nc_bndry,'RDOC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdocwestID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdocwestID,'field','rdoc_west, scalar, series')

    rdonwestID = netcdf.defVar(nc_bndry,'RDON_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdonwestID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdonwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdonwestID,'field','rdon_west, scalar, series')

    rdopwestID = netcdf.defVar(nc_bndry,'RDOP_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdopwestID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rdopwestID,'field','rdop_west, scalar, series')

    redocwestID = netcdf.defVar(nc_bndry,'REDOC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,redocwestID,'long_name','REACTIVE DISSOLVED ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,redocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,redocwestID,'field','redoc_west, scalar, series')

    repocwestID = netcdf.defVar(nc_bndry,'REPOC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,repocwestID,'long_name','REACTIVE PARTICULATE ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,repocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,repocwestID,'field','repoc_west, scalar, series')
    
    rpocwestID = netcdf.defVar(nc_bndry,'RPOC_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpocwestID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpocwestID,'field','rpoc_west, scalar, series')
    
    rponwestID = netcdf.defVar(nc_bndry,'RPON_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rponwestID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponwestID,'field','rpon_west, scalar, series')
    
    rpopwestID = netcdf.defVar(nc_bndry,'RPOP_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpopwestID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopwestID,'field','rpop_west, scalar, series')
    
    salwestID = netcdf.defVar(nc_bndry,'SAL_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,salwestID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,salwestID,'units','psu');
    netcdf.putAtt(nc_bndry,salwestID,'field','sal_west, scalar, series')
    
    sitwestID = netcdf.defVar(nc_bndry,'SIT_west','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,sitwestID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,sitwestID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,sitwestID,'field','sit_west, scalar, series')
end
%EAST
if(dir_flag(4)==1)
    aceastID = netcdf.defVar(nc_bndry,'A_C_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,aceastID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,aceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,aceastID,'field','A_C_east, scalar, series');

    achlceastID = netcdf.defVar(nc_bndry,'A_CHLC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,achlceastID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlceastID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlceastID,'field','A_CHLC_east, scalar, series');

    anceastID = netcdf.defVar(nc_bndry,'A_NC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,anceastID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,anceastID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,anceastID,'field','A_NC_east, scalar, series');

    apceastID = netcdf.defVar(nc_bndry,'A_PC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,apceastID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apceastID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apceastID,'field','A_PC_east, scalar, series');

    bsieastID = netcdf.defVar(nc_bndry,'BSI_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,bsieastID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsieastID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsieastID,'field','bsi_east, scalar, series');

    doeastID = netcdf.defVar(nc_bndry,'DO_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,doeastID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,doeastID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,doeastID,'field','do_east, scalar, series');

    exdoceastID = netcdf.defVar(nc_bndry,'EXDOC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,exdoceastID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdoceastID,'field','exdoc_east, scalar, series');

    ldoceastID = netcdf.defVar(nc_bndry,'LDOC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldoceastID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldoceastID,'field','ldoc_east, scalar, series');

    ldoneastID = netcdf.defVar(nc_bndry,'LDON_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldoneastID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldoneastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldoneastID,'field','ldon_east, scalar, series');

    ldopeastID = netcdf.defVar(nc_bndry,'LDOP_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldopeastID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopeastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopeastID,'field','ldop_east, scalar, series');

    lpoceastID = netcdf.defVar(nc_bndry,'LPOC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpoceastID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpoceastID,'field','lpoc_east, scalar, series');

    lponeastID = netcdf.defVar(nc_bndry,'LPON_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lponeastID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponeastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponeastID,'field','lpon_east, scalar, series');

    lpopeastID = netcdf.defVar(nc_bndry,'LPOP_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpopeastID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopeastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopeastID,'field','lpop_east, scalar, series');

    mavgueastID = netcdf.defVar(nc_bndry,'M_AVGU_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mavgueastID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavgueastID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavgueastID,'field','mavgu_east, scalar, series');

    mceastID = netcdf.defVar(nc_bndry,'M_C_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mceastID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mceastID,'field','mc_east, scalar, series');

    mchlceastID = netcdf.defVar(nc_bndry,'M_CHLC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mchlceastID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlceastID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlceastID,'field','mchlc_east, scalar, series');

    mfceastID = netcdf.defVar(nc_bndry,'M_FC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfceastID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfceastID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfceastID,'field','mfc_east, scalar, series');

    mfchlceastID = netcdf.defVar(nc_bndry,'M_FCHLC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfchlceastID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlceastID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlceastID,'field','mfchlc_east, scalar, series');

    mnceastID = netcdf.defVar(nc_bndry,'M_NC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mnceastID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mnceastID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mnceastID,'field','mnc_east, scalar, series');

    mpceastID = netcdf.defVar(nc_bndry,'M_PC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mpceastID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpceastID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpceastID,'field','mpc_east, scalar, series');

    mfnceastID = netcdf.defVar(nc_bndry,'M_FNC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfnceastID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfnceastID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfnceastID,'field','mfnc_east, scalar, series');

    mfpceastID = netcdf.defVar(nc_bndry,'M_FPC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfpceastID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpceastID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpceastID,'field','mfpc_east, scalar, series');

    nh4teastID = netcdf.defVar(nc_bndry,'NH4T_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,nh4teastID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4teastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4teastID,'field','nh4t_east, scalar, series');

    no23eastID = netcdf.defVar(nc_bndry,'NO23_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,no23eastID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23eastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23eastID,'field','no23_east, scalar, series')

    o2eqeastID = netcdf.defVar(nc_bndry,'O2EQ_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,o2eqeastID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqeastID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqeastID,'field','o2eq_east, scalar, series')

    phyt1eastID = netcdf.defVar(nc_bndry,'PHYT1_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt1eastID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1eastID,'field','phyt1_east, scalar, series')

    phyt2eastID = netcdf.defVar(nc_bndry,'PHYT2_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt2eastID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2eastID,'field','phyt2_east, scalar, series')

    phyt3eastID = netcdf.defVar(nc_bndry,'PHYT3_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt3eastID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3eastID,'field','phyt3_east, scalar, series')

    po4teastID = netcdf.defVar(nc_bndry,'PO4T_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,po4teastID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4teastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4teastID,'field','po4t_east, scalar, series')

    rdoceastID = netcdf.defVar(nc_bndry,'RDOC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdoceastID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdoceastID,'field','rdoc_east, scalar, series')

    rdoneastID = netcdf.defVar(nc_bndry,'RDON_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdoneastID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdoneastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdoneastID,'field','rdon_east, scalar, series')

    rdopeastID = netcdf.defVar(nc_bndry,'RDOP_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdopeastID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopeastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdopeastID,'field','rdop_east, scalar, series')

    redoceastID = netcdf.defVar(nc_bndry,'REDOC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,redoceastID,'long_name','REACTIVE DISSOLVED ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,redoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,redoceastID,'field','redoc_east, scalar, series')

    repoceastID = netcdf.defVar(nc_bndry,'REPOC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,repoceastID,'long_name','REACTIVE PARTICULATE ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,repoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,repoceastID,'field','repoc_east, scalar, series')
    
    rpoceastID = netcdf.defVar(nc_bndry,'RPOC_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpoceastID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpoceastID,'field','rpoc_east, scalar, series')
    
    rponeastID = netcdf.defVar(nc_bndry,'RPON_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rponeastID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponeastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponeastID,'field','rpon_east, scalar, series')
    
    rpopeastID = netcdf.defVar(nc_bndry,'RPOP_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpopeastID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopeastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopeastID,'field','rpop_east, scalar, series')
    
    saleastID = netcdf.defVar(nc_bndry,'SAL_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,saleastID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,saleastID,'units','psu');
    netcdf.putAtt(nc_bndry,saleastID,'field','sal_east, scalar, series')
    
    siteastID = netcdf.defVar(nc_bndry,'SIT_east','float',[yrdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,siteastID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,siteastID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,siteastID,'field','sit_east, scalar, series')
end

%RIVER
if(river>1)
    acriverID = netcdf.defVar(nc_bndry,'river_A_C','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,acriverID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,acriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,acriverID,'field','river_A_C, scalar, series');

    achlcriverID = netcdf.defVar(nc_bndry,'river_A_CHLC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,achlcriverID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlcriverID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlcriverID,'field','river_A_CHLC, scalar, series');

    ancriverID = netcdf.defVar(nc_bndry,'river_A_NC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ancriverID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,ancriverID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,ancriverID,'field','river_A_NC, scalar, series');

    apcriverID = netcdf.defVar(nc_bndry,'river_A_PC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,apcriverID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apcriverID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apcriverID,'field','river_A_PC, scalar, series');

    bsiriverID = netcdf.defVar(nc_bndry,'river_BSI','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,bsiriverID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsiriverID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsiriverID,'field','river_bsi, scalar, series');

    doriverID = netcdf.defVar(nc_bndry,'river_DO','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,doriverID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,doriverID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,doriverID,'field','river_do, scalar, series');

    exdocriverID = netcdf.defVar(nc_bndry,'river_EXDOC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,exdocriverID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdocriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdocriverID,'field','river_exdoc, scalar, series');

    ldocriverID = netcdf.defVar(nc_bndry,'river_LDOC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldocriverID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldocriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldocriverID,'field','river_ldoc, scalar, series');

    ldonriverID = netcdf.defVar(nc_bndry,'river_LDON','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldonriverID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldonriverID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldonriverID,'field','river_ldon, scalar, series');

    ldopriverID = netcdf.defVar(nc_bndry,'river_LDOP','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,ldopriverID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopriverID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopriverID,'field','river_ldop, scalar, series');

    lpocriverID = netcdf.defVar(nc_bndry,'river_LPOC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpocriverID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpocriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpocriverID,'field','river_lpoc, scalar, series');

    lponriverID = netcdf.defVar(nc_bndry,'river_LPON','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lponriverID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponriverID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponriverID,'field','river_lpon, scalar, series');

    lpopriverID = netcdf.defVar(nc_bndry,'river_LPOP','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,lpopriverID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopriverID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopriverID,'field','river_lpop, scalar, series');

    mavguriverID = netcdf.defVar(nc_bndry,'river_M_AVGU','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mavguriverID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavguriverID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavguriverID,'field','river_mavgu, scalar, series');

    mcriverID = netcdf.defVar(nc_bndry,'river_M_C','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mcriverID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mcriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mcriverID,'field','river_mc, scalar, series');

    mchlcriverID = netcdf.defVar(nc_bndry,'river_M_CHLC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mchlcriverID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlcriverID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlcriverID,'field','river_mchlc, scalar, series');

    mfcriverID = netcdf.defVar(nc_bndry,'river_M_FC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfcriverID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfcriverID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfcriverID,'field','river_mfc, scalar, series');

    mfchlcriverID = netcdf.defVar(nc_bndry,'river_M_FCHLC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfchlcriverID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlcriverID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlcriverID,'field','river_mfchlc, scalar, series');

    mncriverID = netcdf.defVar(nc_bndry,'river_M_NC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mncriverID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mncriverID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mncriverID,'field','river_mnc, scalar, series');

    mpcriverID = netcdf.defVar(nc_bndry,'river_M_PC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mpcriverID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpcriverID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpcriverID,'field','river_mpc, scalar, series');

    mfncriverID = netcdf.defVar(nc_bndry,'river_M_FNC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfncriverID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfncriverID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfncriverID,'field','river_mfnc, scalar, series');

    mfpcriverID = netcdf.defVar(nc_bndry,'river_M_FPC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,mfpcriverID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpcriverID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpcriverID,'field','river_mfpc, scalar, series');

    nh4triverID = netcdf.defVar(nc_bndry,'river_NH4T','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,nh4triverID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4triverID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4triverID,'field','river_nh4t, scalar, series');

    no23riverID = netcdf.defVar(nc_bndry,'river_NO23','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,no23riverID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23riverID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23riverID,'field','river_no23, scalar, series')

    o2eqriverID = netcdf.defVar(nc_bndry,'river_O2EQ','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,o2eqriverID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqriverID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqriverID,'field','river_o2eq, scalar, series')

    phyt1riverID = netcdf.defVar(nc_bndry,'river_PHYT1','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt1riverID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1riverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1riverID,'field','river_phyt1, scalar, series')

    phyt2riverID = netcdf.defVar(nc_bndry,'river_PHYT2','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt2riverID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2riverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2riverID,'field','river_phyt2, scalar, series')

    phyt3riverID = netcdf.defVar(nc_bndry,'river_PHYT3','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,phyt3riverID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3riverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3riverID,'field','river_phyt3, scalar, series')

    po4triverID = netcdf.defVar(nc_bndry,'river_PO4T','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,po4triverID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4triverID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4triverID,'field','river_po4t, scalar, series')

    rdocriverID = netcdf.defVar(nc_bndry,'river_RDOC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdocriverID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdocriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdocriverID,'field','river_rdoc, scalar, series')

    rdonriverID = netcdf.defVar(nc_bndry,'river_RDON','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdonriverID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdonriverID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdonriverID,'field','river_rdon, scalar, series')

    rdopriverID = netcdf.defVar(nc_bndry,'river_RDOP','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rdopriverID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopriverID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rdopriverID,'field','river_rdop, scalar, series')

    redocriverID = netcdf.defVar(nc_bndry,'river_REDOC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,redocriverID,'long_name','REACTIVE DISSOLVED ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,redocriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,redocriverID,'field','river_redoc, scalar, series')

    repocriverID = netcdf.defVar(nc_bndry,'river_REPOC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,repocriverID,'long_name','REACTIVE PARTICULATE ORGANIC CARBON - CSO/WWTP');
    netcdf.putAtt(nc_bndry,repocriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,repocriverID,'field','river_repoc, scalar, series')
    
    rpocriverID = netcdf.defVar(nc_bndry,'river_RPOC','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpocriverID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpocriverID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpocriverID,'field','river_rpoc, scalar, series')
    
    rponriverID = netcdf.defVar(nc_bndry,'river_RPON','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rponriverID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponriverID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponriverID,'field','river_rpon, scalar, series')
    
    rpopriverID = netcdf.defVar(nc_bndry,'river_RPOP','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,rpopriverID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopriverID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopriverID,'field','river_rpop, scalar, series')
    
    salriverID = netcdf.defVar(nc_bndry,'river_SAL','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,salriverID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,salriverID,'units','psu');
    netcdf.putAtt(nc_bndry,salriverID,'field','river_sal, scalar, series')
    
    sitriverID = netcdf.defVar(nc_bndry,'river_SIT','float',[rdimID srdimID tdimID]);
    netcdf.putAtt(nc_bndry,sitriverID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,sitriverID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,sitriverID,'field','river_sit, scalar, series')
end
netcdf.close(nc_bndry)

end

