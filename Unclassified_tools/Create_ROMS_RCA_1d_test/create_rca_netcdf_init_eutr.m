function create_rca_netcdf_init_eutr(fn,gn,t)
% Generate a blank netcdf file containing RCA EUTRO initial conditions.
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
nc_init=netcdf.create(fn,'clobber');
if isempty(nc_init), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by cyr on ' datestr(now)]);
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'type', 'Nutrients data from Florida Water Atlas');
%% Dimensions:

disp(' ## Defining Dimensions...')
 
xrdimID = netcdf.defDim(nc_init,'x_r',x_r);
yrdimID = netcdf.defDim(nc_init,'y_r',y_r);
srdimID = netcdf.defDim(nc_init,'s_r',s);
tdimID = netcdf.defDim(nc_init,'time',t);

%% Variables and attributes:
tID = netcdf.defVar(nc_init,'TIME','float',tdimID);
netcdf.putAtt(nc_init,tID,'long_name','initial condition time');
netcdf.putAtt(nc_init,tID,'units','days');
netcdf.putAtt(nc_init,tID,'field','TIME, scalar, series');

acID = netcdf.defVar(nc_init,'A_C','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,acID,'long_name','MX model prey');
netcdf.putAtt(nc_init,acID,'units','mg C L-1');
netcdf.putAtt(nc_init,acID,'field','A_C, scalar, series');

achlcID = netcdf.defVar(nc_init,'A_CHLC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,achlcID,'long_name','MX model CHLC/C of Prey');
netcdf.putAtt(nc_init,achlcID,'units','gChl/gC');
netcdf.putAtt(nc_init,achlcID,'field','A_CHLC, scalar, series');

ancID = netcdf.defVar(nc_init,'A_NC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,ancID,'long_name','MX model N:C of Prey');
netcdf.putAtt(nc_init,ancID,'units','gN/gC');
netcdf.putAtt(nc_init,ancID,'field','A_NC, scalar, series');

apcID = netcdf.defVar(nc_init,'A_PC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,apcID,'long_name','MX model P:C of Prey');
netcdf.putAtt(nc_init,apcID,'units','gP/gC');
netcdf.putAtt(nc_init,apcID,'field','A_PC, scalar, series');

bsiID = netcdf.defVar(nc_init,'BSI','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,bsiID,'long_name','BIOGENIC SILICA');
netcdf.putAtt(nc_init,bsiID,'units','mg SI L-1');
netcdf.putAtt(nc_init,bsiID,'field','bsi, scalar, series');

doID = netcdf.defVar(nc_init,'DO','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,doID,'long_name','DISSOLVED OXYGEN');
netcdf.putAtt(nc_init,doID,'units','mg O2 L-1');
netcdf.putAtt(nc_init,doID,'field','do, scalar, series');

exdocID = netcdf.defVar(nc_init,'EXDOC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,exdocID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
netcdf.putAtt(nc_init,exdocID,'units','mg C L-1');
netcdf.putAtt(nc_init,exdocID,'field','exdoc, scalar, series');

ldocID = netcdf.defVar(nc_init,'LDOC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,ldocID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
netcdf.putAtt(nc_init,ldocID,'units','mg C L-1');
netcdf.putAtt(nc_init,ldocID,'field','ldoc, scalar, series');

ldonID = netcdf.defVar(nc_init,'LDON','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,ldonID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
netcdf.putAtt(nc_init,ldonID,'units','mg N L-1');
netcdf.putAtt(nc_init,ldonID,'field','ldon, scalar, series');

ldopID = netcdf.defVar(nc_init,'LDOP','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,ldopID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
netcdf.putAtt(nc_init,ldopID,'units','mg P L-1');
netcdf.putAtt(nc_init,ldopID,'field','ldop, scalar, series');


lpocID = netcdf.defVar(nc_init,'LPOC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,lpocID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
netcdf.putAtt(nc_init,lpocID,'units','mg C L-1');
netcdf.putAtt(nc_init,lpocID,'field','lpoc, scalar, series');

lponID = netcdf.defVar(nc_init,'LPON','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,lponID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
netcdf.putAtt(nc_init,lponID,'units','mg N L-1');
netcdf.putAtt(nc_init,lponID,'field','lpon, scalar, series');

lpopID = netcdf.defVar(nc_init,'LPOP','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,lpopID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
netcdf.putAtt(nc_init,lpopID,'units','mg P L-1');
netcdf.putAtt(nc_init,lpopID,'field','lpop, scalar, series');

mavguID = netcdf.defVar(nc_init,'M_AVGU','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mavguID,'long_name','MX model mixotroph average growth rate');
netcdf.putAtt(nc_init,mavguID,'units','C/C/d');
netcdf.putAtt(nc_init,mavguID,'field','mavgu, scalar, series');

mcID = netcdf.defVar(nc_init,'M_C','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mcID,'long_name','MX model core mixotroph body biomass carbon');
netcdf.putAtt(nc_init,mcID,'units','mg C L-1');
netcdf.putAtt(nc_init,mcID,'field','mc, scalar, series');

mchlcID = netcdf.defVar(nc_init,'M_CHLC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mchlcID,'long_name','MX model mixotroph core Chla:C');
netcdf.putAtt(nc_init,mchlcID,'units','gChl/gC');
netcdf.putAtt(nc_init,mchlcID,'field','mchlc, scalar, series');

mfcID = netcdf.defVar(nc_init,'M_FC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mfcID,'long_name','MX model mixotroph gut relative to body biomass');
netcdf.putAtt(nc_init,mfcID,'units','gC/gC');
netcdf.putAtt(nc_init,mfcID,'field','mfc, scalar, series');

mfchlcID = netcdf.defVar(nc_init,'M_FCHLC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mfchlcID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
netcdf.putAtt(nc_init,mfchlcID,'units','gChl in the gut /gC');
netcdf.putAtt(nc_init,mfchlcID,'field','mfchlc, scalar, series');

mncID = netcdf.defVar(nc_init,'M_NC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mncID,'long_name','MX model mixotroph N:C');
netcdf.putAtt(nc_init,mncID,'units','gN/gC');
netcdf.putAtt(nc_init,mncID,'field','mnc, scalar, series');

mpcID = netcdf.defVar(nc_init,'M_PC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mpcID,'long_name','MX model mixotroph P:C');
netcdf.putAtt(nc_init,mpcID,'units','gP/gC');
netcdf.putAtt(nc_init,mpcID,'field','mpc, scalar, series');

mfncID = netcdf.defVar(nc_init,'M_FNC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mfncID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
netcdf.putAtt(nc_init,mfncID,'units','gN/gC');
netcdf.putAtt(nc_init,mfncID,'field','mfnc, scalar, series');

mfpcID = netcdf.defVar(nc_init,'M_FPC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,mfpcID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
netcdf.putAtt(nc_init,mfpcID,'units','gP/gC');
netcdf.putAtt(nc_init,mfpcID,'field','mfpc, scalar, series');

nh4tID = netcdf.defVar(nc_init,'NH4T','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,nh4tID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
netcdf.putAtt(nc_init,nh4tID,'units','mg N L-1');
netcdf.putAtt(nc_init,nh4tID,'field','nh4t, scalar, series');

no23ID = netcdf.defVar(nc_init,'NO23','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,no23ID,'long_name','NITRITE + NITRATE');
netcdf.putAtt(nc_init,no23ID,'units','mg N L-1');
netcdf.putAtt(nc_init,no23ID,'field','no23, scalar, series')

o2eqID = netcdf.defVar(nc_init,'O2EQ','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,o2eqID,'long_name','AQUEOUS SOD');
netcdf.putAtt(nc_init,o2eqID,'units','mg O2 L-1');
netcdf.putAtt(nc_init,o2eqID,'field','o2eq, scalar, series')

phyt1ID = netcdf.defVar(nc_init,'PHYT1','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,phyt1ID,'long_name','WINTER DIATOMS (PHYT1)');
netcdf.putAtt(nc_init,phyt1ID,'units','mg C L-1');
netcdf.putAtt(nc_init,phyt1ID,'field','phyt1, scalar, series')

phyt2ID = netcdf.defVar(nc_init,'PHYT2','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,phyt2ID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
netcdf.putAtt(nc_init,phyt2ID,'units','mg C L-1');
netcdf.putAtt(nc_init,phyt2ID,'field','phyt2, scalar, series')

phyt3ID = netcdf.defVar(nc_init,'PHYT3','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,phyt3ID,'long_name','FALL ASSEMBLAGE (PHYT3)');
netcdf.putAtt(nc_init,phyt3ID,'units','mg C L-1');
netcdf.putAtt(nc_init,phyt3ID,'field','phyt3, scalar, series')

po4tID = netcdf.defVar(nc_init,'PO4T','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,po4tID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
netcdf.putAtt(nc_init,po4tID,'units','mg P L-1');
netcdf.putAtt(nc_init,po4tID,'field','po4t, scalar, series')

rdocID = netcdf.defVar(nc_init,'RDOC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,rdocID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
netcdf.putAtt(nc_init,rdocID,'units','mg C L-1');
netcdf.putAtt(nc_init,rdocID,'field','rdoc, scalar, series')

rdonID = netcdf.defVar(nc_init,'RDON','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,rdonID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
netcdf.putAtt(nc_init,rdonID,'units','mg N L-1');
netcdf.putAtt(nc_init,rdonID,'field','rdon, scalar, series')

rdopID = netcdf.defVar(nc_init,'RDOP','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,rdopID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
netcdf.putAtt(nc_init,rdopID,'units','mg P L-1');
netcdf.putAtt(nc_init,rdopID,'field','rdop, scalar, series')

redocID = netcdf.defVar(nc_init,'REDOC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,redocID,'long_name','REACTIVE DISSOLVED ORGANIC CARBON - CSO/WWTP');
netcdf.putAtt(nc_init,redocID,'units','mg C L-1');
netcdf.putAtt(nc_init,redocID,'field','redoc, scalar, series')

repocID = netcdf.defVar(nc_init,'REPOC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,repocID,'long_name','REACTIVE PARTICULATE ORGANIC CARBON - CSO/WWTP');
netcdf.putAtt(nc_init,repocID,'units','mg C L-1');
netcdf.putAtt(nc_init,repocID,'field','repoc, scalar, series')

rpocID = netcdf.defVar(nc_init,'RPOC','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,rpocID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
netcdf.putAtt(nc_init,rpocID,'units','mg C L-1');
netcdf.putAtt(nc_init,rpocID,'field','rpoc, scalar, series')

rponID = netcdf.defVar(nc_init,'RPON','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,rponID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
netcdf.putAtt(nc_init,rponID,'units','mg N L-1');
netcdf.putAtt(nc_init,rponID,'field','rpon, scalar, series')

rpopID = netcdf.defVar(nc_init,'RPOP','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,rpopID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
netcdf.putAtt(nc_init,rpopID,'units','mg P L-1');
netcdf.putAtt(nc_init,rpopID,'field','rpop, scalar, series')

salID = netcdf.defVar(nc_init,'SAL','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,salID,'long_name','SALINITY');
netcdf.putAtt(nc_init,salID,'units','psu');
netcdf.putAtt(nc_init,salID,'field','sal, scalar, series')

sitID = netcdf.defVar(nc_init,'SIT','float',[xrdimID yrdimID srdimID tdimID]);
netcdf.putAtt(nc_init,sitID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
netcdf.putAtt(nc_init,sitID,'units','mg SI L-1');
netcdf.putAtt(nc_init,sitID,'field','sit, scalar, series')

netcdf.close(nc_init)

