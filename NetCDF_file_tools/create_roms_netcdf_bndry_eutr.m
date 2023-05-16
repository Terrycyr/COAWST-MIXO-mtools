function create_roms_netcdf_bndry_eutr(fn,gn,t,dir_flag)
% Generate a blank netcdf file containing RCA EUTRO boundary conditions.
% Yuren Chen, 2021/06/21, Guangzhou
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
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by cyr on ' datestr(now)]);
netcdf.putAtt(nc_bndry,netcdf.getConstant('NC_GLOBAL'),'type', 'Nutrients data from Florida Water Atlas');
%% Dimensions:

disp(' ## Defining Dimensions...')
 
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

biotdimID = netcdf.defDim(nc_bndry,'bio_time',t);
SALtdimID = netcdf.defDim(nc_bndry,'SAL_time',t);
PHYT1tdimID = netcdf.defDim(nc_bndry,'PHYT1_time',t);
PHYT2tdimID = netcdf.defDim(nc_bndry,'PHYT2_time',t);
PHYT3tdimID = netcdf.defDim(nc_bndry,'PHYT3_time',t);
RPOPtdimID = netcdf.defDim(nc_bndry,'RPOP_time',t);
LPOPtdimID = netcdf.defDim(nc_bndry,'LPOP_time',t);
RDOPtdimID = netcdf.defDim(nc_bndry,'RDOP_time',t);
LDOPtdimID = netcdf.defDim(nc_bndry,'LDOP_time',t);
PO4TtdimID = netcdf.defDim(nc_bndry,'PO4T_time',t);
RPONtdimID = netcdf.defDim(nc_bndry,'RPON_time',t);
LPONtdimID = netcdf.defDim(nc_bndry,'LPON_time',t);
RDONtdimID = netcdf.defDim(nc_bndry,'RDON_time',t);
LDONtdimID = netcdf.defDim(nc_bndry,'LDON_time',t);
NH4TtdimID = netcdf.defDim(nc_bndry,'NH4T_time',t);
NO23tdimID = netcdf.defDim(nc_bndry,'NO23_time',t);
BSItdimID = netcdf.defDim(nc_bndry,'BSI_time',t);
SITtdimID = netcdf.defDim(nc_bndry,'SIT_time',t);
RPOCtdimID = netcdf.defDim(nc_bndry,'RPOC_time',t);
LPOCtdimID = netcdf.defDim(nc_bndry,'LPOC_time',t);
RDOCtdimID = netcdf.defDim(nc_bndry,'RDOC_time',t);
LDOCtdimID = netcdf.defDim(nc_bndry,'LDOC_time',t);
EXDOCtdimID = netcdf.defDim(nc_bndry,'EXDOC_time',t);
DICtdimID = netcdf.defDim(nc_bndry,'DIC_time',t);
CACO3tdimID = netcdf.defDim(nc_bndry,'CACO3_time',t);
TAtdimID = netcdf.defDim(nc_bndry,'TA_time',t);
O2EQtdimID = netcdf.defDim(nc_bndry,'O2EQ_time',t);
DOtdimID = netcdf.defDim(nc_bndry,'DO_time',t);

A_CtdimID = netcdf.defDim(nc_bndry,'A_C_time',t);
A_CHLCtdimID = netcdf.defDim(nc_bndry,'A_CHLC_time',t);
A_NCtdimID = netcdf.defDim(nc_bndry,'A_NC_time',t);
A_PCtdimID = netcdf.defDim(nc_bndry,'A_PC_time',t);
M_AVGCUtdimID = netcdf.defDim(nc_bndry,'M_AVGCU_time',t);
M_CtdimID = netcdf.defDim(nc_bndry,'M_C_time',t);
M_CHLCtdimID = netcdf.defDim(nc_bndry,'M_CHLC_time',t);
M_FCtdimID = netcdf.defDim(nc_bndry,'M_FC_time',t);
M_FCHLCtdimID = netcdf.defDim(nc_bndry,'M_FCHLC_time',t);
M_NCtdimID = netcdf.defDim(nc_bndry,'M_NC_time',t);
M_PCtdimID = netcdf.defDim(nc_bndry,'M_PC_time',t);
MFNCtdimID = netcdf.defDim(nc_bndry,'MFNC_time',t);
MFPCtdimID = netcdf.defDim(nc_bndry,'MFPC_time',t);

vname = {'bio','SAL','PHYT1','PHYT2','PHYT3','RPOP','LPOP','RDOP','LDOP','PO4T'...
    ,'RPON','LPON','RDON','LDON','NH4T','NO23','BSI','SIT','RPOC','LPOC'...
    ,'RDOC','LDOC','EXDOC','O2EQ','DO','A_C','A_CHLC','A_NC','A_PC','M_AVGCU'...
    ,'M_C','M_CHLC','M_FC','M_FCHLC','M_NC','M_PC','MFNC','MFPC','DIC','TA','CACO3','ZOO1','ZOO2'};

%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
%
for ivar = 1:40
    eval(['tID = netcdf.defVar(nc_bndry,''',vname{ivar},'_time'',''float'',',vname{ivar},'tdimID);']);
    eval(['netcdf.putAtt(nc_bndry,tID,''long_name'',''',vname{ivar},'_time'');']);
    netcdf.putAtt(nc_bndry,tID,'units','days');
    eval(['netcdf.putAtt(nc_bndry,tID,''field'',''',vname{ivar},'_time, scalar, series'');']);
end


%South
if(dir_flag(1)==1)
    acsouthID = netcdf.defVar(nc_bndry,'A_C_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,acsouthID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,acsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,acsouthID,'field','A_C, scalar, series');

    achlcsouthID = netcdf.defVar(nc_bndry,'A_CHLC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,achlcsouthID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlcsouthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlcsouthID,'field','A_CHLC, scalar, series');

    ancsouthID = netcdf.defVar(nc_bndry,'A_NC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ancsouthID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,ancsouthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,ancsouthID,'field','A_NC, scalar, series');

    apcsouthID = netcdf.defVar(nc_bndry,'A_PC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,apcsouthID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apcsouthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apcsouthID,'field','A_PC, scalar, series');

    bsisouthID = netcdf.defVar(nc_bndry,'BSI_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,bsisouthID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsisouthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsisouthID,'field','BSI, scalar, series');

    dosouthID = netcdf.defVar(nc_bndry,'DO_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,dosouthID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,dosouthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,dosouthID,'field','DO, scalar, series');

    exdocsouthID = netcdf.defVar(nc_bndry,'EXDOC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,exdocsouthID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdocsouthID,'field','EXDOC, scalar, series');

    ldocsouthID = netcdf.defVar(nc_bndry,'LDOC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldocsouthID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldocsouthID,'field','LDOC, scalar, series');

    ldonsouthID = netcdf.defVar(nc_bndry,'LDON_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldonsouthID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldonsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldonsouthID,'field','LDON, scalar, series');

    ldopsouthID = netcdf.defVar(nc_bndry,'LDOP_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldopsouthID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopsouthID,'field','LDOP, scalar, series');


    lpocsouthID = netcdf.defVar(nc_bndry,'LPOC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpocsouthID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpocsouthID,'field','LPOC, scalar, series');

    lponsouthID = netcdf.defVar(nc_bndry,'LPON_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lponsouthID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponsouthID,'field','LPON, scalar, series');

    lpopsouthID = netcdf.defVar(nc_bndry,'LPOP_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpopsouthID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopsouthID,'field','LPOP, scalar, series');

    mavgusouthID = netcdf.defVar(nc_bndry,'M_AVGCU_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mavgusouthID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavgusouthID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavgusouthID,'field','MAVGCU, scalar, series');

    mcsouthID = netcdf.defVar(nc_bndry,'M_C_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mcsouthID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mcsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mcsouthID,'field','MC, scalar, series');

    mchlcsouthID = netcdf.defVar(nc_bndry,'M_CHLC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mchlcsouthID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlcsouthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlcsouthID,'field','MCHLC, scalar, series');

    mfcsouthID = netcdf.defVar(nc_bndry,'M_FC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfcsouthID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfcsouthID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfcsouthID,'field','MFC, scalar, series');

    mfchlcsouthID = netcdf.defVar(nc_bndry,'M_FCHLC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfchlcsouthID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlcsouthID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlcsouthID,'field','MFCHLC, scalar, series');

    mncsouthID = netcdf.defVar(nc_bndry,'M_NC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mncsouthID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mncsouthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mncsouthID,'field','MNC, scalar, series');

    mpcsouthID = netcdf.defVar(nc_bndry,'M_PC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mpcsouthID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpcsouthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpcsouthID,'field','MPC, scalar, series');

    mfncsouthID = netcdf.defVar(nc_bndry,'MFNC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfncsouthID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfncsouthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfncsouthID,'field','MFNC, scalar, series');

    mfpcsouthID = netcdf.defVar(nc_bndry,'MFPC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfpcsouthID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpcsouthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpcsouthID,'field','MFPC, scalar, series');

    nh4tsouthID = netcdf.defVar(nc_bndry,'NH4T_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,nh4tsouthID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4tsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4tsouthID,'field','NH4T, scalar, series');

    no23southID = netcdf.defVar(nc_bndry,'NO23_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,no23southID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23southID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23southID,'field','NO23, scalar, series')

    o2eqsouthID = netcdf.defVar(nc_bndry,'O2EQ_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,o2eqsouthID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqsouthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqsouthID,'field','O2EQ, scalar, series')

    phyt1southID = netcdf.defVar(nc_bndry,'PHYT1_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt1southID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1southID,'field','PHYT1, scalar, series')

    phyt2southID = netcdf.defVar(nc_bndry,'PHYT2_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt2southID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2southID,'field','PHYT2, scalar, series')

    phyt3southID = netcdf.defVar(nc_bndry,'PHYT3_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt3southID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3southID,'field','PHYT3, scalar, series')

    zoo1southID = netcdf.defVar(nc_bndry,'ZOO1_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo1southID,'long_name','SMALL ZOOPLANKTON (ZOO1)');
    netcdf.putAtt(nc_bndry,zoo1southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo1southID,'field','ZOO1, scalar, series')

    zoo2southID = netcdf.defVar(nc_bndry,'ZOO2_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo2southID,'long_name','LARGE ZOOPLANKTON (ZOO2)');
    netcdf.putAtt(nc_bndry,zoo2southID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo2southID,'field','ZOO2, scalar, series')

    po4tsouthID = netcdf.defVar(nc_bndry,'PO4T_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,po4tsouthID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4tsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4tsouthID,'field','PO4T, scalar, series')

    rdocsouthID = netcdf.defVar(nc_bndry,'RDOC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdocsouthID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdocsouthID,'field','RDOC, scalar, series')

    rdonsouthID = netcdf.defVar(nc_bndry,'RDON_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdonsouthID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdonsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdonsouthID,'field','RDON, scalar, series')

    rdopsouthID = netcdf.defVar(nc_bndry,'RDOP_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdopsouthID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rdopsouthID,'field','RDOP, scalar, series')

    TAsouthID = netcdf.defVar(nc_bndry,'TA_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,TAsouthID,'long_name','TOTAL ALKALINITY');
    netcdf.putAtt(nc_bndry,TAsouthID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,TAsouthID,'field','TA, scalar, series')

    DICsouthID = netcdf.defVar(nc_bndry,'DIC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,DICsouthID,'long_name','DISSOLVED INORGANIC CARBON');
    netcdf.putAtt(nc_bndry,DICsouthID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,DICsouthID,'field','DIC, scalar, series')

    CACO3southID = netcdf.defVar(nc_bndry,'CACO3_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,CACO3southID,'long_name','CALCIUM CARBONATE');
    netcdf.putAtt(nc_bndry,CACO3southID,'units','mmol kg-1');
    netcdf.putAtt(nc_bndry,CACO3southID,'field','CACO3, scalar, series')
    
    rpocsouthID = netcdf.defVar(nc_bndry,'RPOC_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpocsouthID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpocsouthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpocsouthID,'field','RPOC, scalar, series')
    
    rponsouthID = netcdf.defVar(nc_bndry,'RPON_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rponsouthID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponsouthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponsouthID,'field','RPON, scalar, series')
    
    rpopsouthID = netcdf.defVar(nc_bndry,'RPOP_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpopsouthID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopsouthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopsouthID,'field','RPOP, scalar, series')
    
    salsouthID = netcdf.defVar(nc_bndry,'SAL_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,salsouthID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,salsouthID,'units','psu');
    netcdf.putAtt(nc_bndry,salsouthID,'field','SAL, scalar, series')
    
    sitsouthID = netcdf.defVar(nc_bndry,'SIT_south','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,sitsouthID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,sitsouthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,sitsouthID,'field','SIT, scalar, series')
end

%NORTH
if(dir_flag(2)==1)
    acnorthID = netcdf.defVar(nc_bndry,'A_C_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,acnorthID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,acnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,acnorthID,'field','A_C, scalar, series');

    achlcnorthID = netcdf.defVar(nc_bndry,'A_CHLC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,achlcnorthID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlcnorthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlcnorthID,'field','A_CHLC, scalar, series');

    ancnorthID = netcdf.defVar(nc_bndry,'A_NC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ancnorthID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,ancnorthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,ancnorthID,'field','A_NC, scalar, series');

    apcnorthID = netcdf.defVar(nc_bndry,'A_PC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,apcnorthID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apcnorthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apcnorthID,'field','A_PC, scalar, series');

    bsinorthID = netcdf.defVar(nc_bndry,'BSI_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,bsinorthID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsinorthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsinorthID,'field','BSI, scalar, series');

    donorthID = netcdf.defVar(nc_bndry,'DO_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,donorthID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,donorthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,donorthID,'field','DO, scalar, series');

    exdocnorthID = netcdf.defVar(nc_bndry,'EXDOC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,exdocnorthID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdocnorthID,'field','EXDOC, scalar, series');

    ldocnorthID = netcdf.defVar(nc_bndry,'LDOC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldocnorthID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldocnorthID,'field','LDOC, scalar, series');

    ldonnorthID = netcdf.defVar(nc_bndry,'LDON_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldonnorthID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldonnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldonnorthID,'field','LDON, scalar, series');

    ldopnorthID = netcdf.defVar(nc_bndry,'LDOP_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldopnorthID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopnorthID,'field','LDOP, scalar, series');

    lpocnorthID = netcdf.defVar(nc_bndry,'LPOC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpocnorthID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpocnorthID,'field','LPOC, scalar, series');

    lponnorthID = netcdf.defVar(nc_bndry,'LPON_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lponnorthID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponnorthID,'field','LPON, scalar, series');

    lpopnorthID = netcdf.defVar(nc_bndry,'LPOP_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpopnorthID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopnorthID,'field','LPOP, scalar, series');

    mavgunorthID = netcdf.defVar(nc_bndry,'M_AVGCU_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mavgunorthID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavgunorthID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavgunorthID,'field','MAVGCU, scalar, series');

    mcnorthID = netcdf.defVar(nc_bndry,'M_C_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mcnorthID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mcnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mcnorthID,'field','MC, scalar, series');

    mchlcnorthID = netcdf.defVar(nc_bndry,'M_CHLC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mchlcnorthID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlcnorthID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlcnorthID,'field','MCHLC, scalar, series');

    mfcnorthID = netcdf.defVar(nc_bndry,'M_FC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfcnorthID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfcnorthID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfcnorthID,'field','MFC, scalar, series');

    mfchlcnorthID = netcdf.defVar(nc_bndry,'M_FCHLC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfchlcnorthID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlcnorthID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlcnorthID,'field','MFCHLC, scalar, series');

    mncnorthID = netcdf.defVar(nc_bndry,'M_NC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mncnorthID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mncnorthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mncnorthID,'field','MNC, scalar, series');

    mpcnorthID = netcdf.defVar(nc_bndry,'M_PC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mpcnorthID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpcnorthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpcnorthID,'field','MPC, scalar, series');

    mfncnorthID = netcdf.defVar(nc_bndry,'MFNC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfncnorthID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfncnorthID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfncnorthID,'field','MFNC, scalar, series');

    mfpcnorthID = netcdf.defVar(nc_bndry,'MFPC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfpcnorthID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpcnorthID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpcnorthID,'field','MFPC, scalar, series');

    nh4tnorthID = netcdf.defVar(nc_bndry,'NH4T_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,nh4tnorthID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4tnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4tnorthID,'field','NH4T, scalar, series');

    no23northID = netcdf.defVar(nc_bndry,'NO23_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,no23northID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23northID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23northID,'field','NO23, scalar, series')

    o2eqnorthID = netcdf.defVar(nc_bndry,'O2EQ_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,o2eqnorthID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqnorthID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqnorthID,'field','O2EQ, scalar, series')

    phyt1northID = netcdf.defVar(nc_bndry,'PHYT1_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt1northID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1northID,'field','PHYT1, scalar, series')

    phyt2northID = netcdf.defVar(nc_bndry,'PHYT2_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt2northID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2northID,'field','PHYT2, scalar, series')

    phyt3northID = netcdf.defVar(nc_bndry,'PHYT3_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt3northID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3northID,'field','PHYT3, scalar, series')

    zoo1northID = netcdf.defVar(nc_bndry,'ZOO1_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo1northID,'long_name','SMALL ZOOPLANKTON (ZOO1)');
    netcdf.putAtt(nc_bndry,zoo1northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo1northID,'field','ZOO1, scalar, series')

    zoo2northID = netcdf.defVar(nc_bndry,'ZOO2_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo2northID,'long_name','LARGE ZOOPLANKTON (ZOO2)');
    netcdf.putAtt(nc_bndry,zoo2northID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo2northID,'field','ZOO2, scalar, series')

    po4tnorthID = netcdf.defVar(nc_bndry,'PO4T_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,po4tnorthID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4tnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4tnorthID,'field','PO4T, scalar, series')

    rdocnorthID = netcdf.defVar(nc_bndry,'RDOC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdocnorthID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdocnorthID,'field','RDOC, scalar, series')

    rdonnorthID = netcdf.defVar(nc_bndry,'RDON_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdonnorthID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdonnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdonnorthID,'field','RDON, scalar, series')

    rdopnorthID = netcdf.defVar(nc_bndry,'RDOP_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdopnorthID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rdopnorthID,'field','RDOP, scalar, series')

    TAnorthID = netcdf.defVar(nc_bndry,'TA_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,TAnorthID,'long_name','TOTAL ALKALINITY');
    netcdf.putAtt(nc_bndry,TAnorthID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,TAnorthID,'field','TA, scalar, series')

    DICnorthID = netcdf.defVar(nc_bndry,'DIC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,DICnorthID,'long_name','DISSOLVED INORGANIC CARBON');
    netcdf.putAtt(nc_bndry,DICnorthID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,DICnorthID,'field','DIC, scalar, series')

    CACO3northID = netcdf.defVar(nc_bndry,'CACO3_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,CACO3northID,'long_name','CALCIUM CARBONATE');
    netcdf.putAtt(nc_bndry,CACO3northID,'units','mmol kg-1');
    netcdf.putAtt(nc_bndry,CACO3northID,'field','CACO3, scalar, series')
    
    rpocnorthID = netcdf.defVar(nc_bndry,'RPOC_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpocnorthID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpocnorthID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpocnorthID,'field','RPOC, scalar, series')
    
    rponnorthID = netcdf.defVar(nc_bndry,'RPON_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rponnorthID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponnorthID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponnorthID,'field','RPON, scalar, series')
    
    rpopnorthID = netcdf.defVar(nc_bndry,'RPOP_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpopnorthID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopnorthID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopnorthID,'field','RPOP, scalar, series')
    
    salnorthID = netcdf.defVar(nc_bndry,'SAL_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,salnorthID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,salnorthID,'units','psu');
    netcdf.putAtt(nc_bndry,salnorthID,'field','SAL, scalar, series')
    
    sitnorthID = netcdf.defVar(nc_bndry,'SIT_north','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,sitnorthID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,sitnorthID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,sitnorthID,'field','SIT, scalar, series')
end
%WEST
if(dir_flag(3)==1)
    acwestID = netcdf.defVar(nc_bndry,'A_C_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,acwestID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,acwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,acwestID,'field','A_C, scalar, series');

    achlcwestID = netcdf.defVar(nc_bndry,'A_CHLC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,achlcwestID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlcwestID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlcwestID,'field','A_CHLC, scalar, series');

    ancwestID = netcdf.defVar(nc_bndry,'A_NC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ancwestID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,ancwestID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,ancwestID,'field','A_NC, scalar, series');

    apcwestID = netcdf.defVar(nc_bndry,'A_PC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,apcwestID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apcwestID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apcwestID,'field','A_PC, scalar, series');

    bsiwestID = netcdf.defVar(nc_bndry,'BSI_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,bsiwestID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsiwestID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsiwestID,'field','BSI, scalar, series');

    dowestID = netcdf.defVar(nc_bndry,'DO_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,dowestID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,dowestID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,dowestID,'field','DO, scalar, series');

    exdocwestID = netcdf.defVar(nc_bndry,'EXDOC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,exdocwestID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdocwestID,'field','EXDOC, scalar, series');

    ldocwestID = netcdf.defVar(nc_bndry,'LDOC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldocwestID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldocwestID,'field','LDOC, scalar, series');

    ldonwestID = netcdf.defVar(nc_bndry,'LDON_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldonwestID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldonwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldonwestID,'field','LDON, scalar, series');

    ldopwestID = netcdf.defVar(nc_bndry,'LDOP_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldopwestID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopwestID,'field','LDOP, scalar, series');

    lpocwestID = netcdf.defVar(nc_bndry,'LPOC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpocwestID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpocwestID,'field','LPOC, scalar, series');

    lponwestID = netcdf.defVar(nc_bndry,'LPON_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lponwestID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponwestID,'field','LPON, scalar, series');

    lpopwestID = netcdf.defVar(nc_bndry,'LPOP_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpopwestID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopwestID,'field','LPOP, scalar, series');

    mavguwestID = netcdf.defVar(nc_bndry,'M_AVGCU_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mavguwestID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavguwestID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavguwestID,'field','MAVGCU, scalar, series');

    mcwestID = netcdf.defVar(nc_bndry,'M_C_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mcwestID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mcwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mcwestID,'field','MC, scalar, series');

    mchlcwestID = netcdf.defVar(nc_bndry,'M_CHLC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mchlcwestID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlcwestID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlcwestID,'field','MCHLC, scalar, series');

    mfcwestID = netcdf.defVar(nc_bndry,'M_FC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfcwestID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfcwestID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfcwestID,'field','MFC, scalar, series');

    mfchlcwestID = netcdf.defVar(nc_bndry,'M_FCHLC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfchlcwestID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlcwestID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlcwestID,'field','MFCHLC, scalar, series');

    mncwestID = netcdf.defVar(nc_bndry,'M_NC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mncwestID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mncwestID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mncwestID,'field','MNC, scalar, series');

    mpcwestID = netcdf.defVar(nc_bndry,'M_PC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mpcwestID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpcwestID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpcwestID,'field','MPC, scalar, series');

    mfncwestID = netcdf.defVar(nc_bndry,'MFNC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfncwestID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfncwestID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfncwestID,'field','MFNC, scalar, series');

    mfpcwestID = netcdf.defVar(nc_bndry,'MFPC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfpcwestID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpcwestID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpcwestID,'field','MFPC, scalar, series');

    nh4twestID = netcdf.defVar(nc_bndry,'NH4T_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,nh4twestID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4twestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4twestID,'field','NH4T, scalar, series');

    no23westID = netcdf.defVar(nc_bndry,'NO23_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,no23westID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23westID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23westID,'field','NO23, scalar, series')

    o2eqwestID = netcdf.defVar(nc_bndry,'O2EQ_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,o2eqwestID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqwestID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqwestID,'field','O2EQ, scalar, series')

    phyt1westID = netcdf.defVar(nc_bndry,'PHYT1_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt1westID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1westID,'field','PHYT1, scalar, series')

    phyt2westID = netcdf.defVar(nc_bndry,'PHYT2_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt2westID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2westID,'field','PHYT2, scalar, series')

    phyt3westID = netcdf.defVar(nc_bndry,'PHYT3_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt3westID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3westID,'field','PHYT3, scalar, series')

    zoo1westID = netcdf.defVar(nc_bndry,'ZOO1_west','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo1westID,'long_name','SMALL ZOOPLANKTON (ZOO1)');
    netcdf.putAtt(nc_bndry,zoo1westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo1westID,'field','ZOO1, scalar, series')

    zoo2westID = netcdf.defVar(nc_bndry,'ZOO2_west','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo2westID,'long_name','LARGE ZOOPLANKTON (ZOO2)');
    netcdf.putAtt(nc_bndry,zoo2westID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo2westID,'field','ZOO2, scalar, series')

    po4twestID = netcdf.defVar(nc_bndry,'PO4T_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,po4twestID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4twestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4twestID,'field','PO4T, scalar, series')

    rdocwestID = netcdf.defVar(nc_bndry,'RDOC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdocwestID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdocwestID,'field','RDOC, scalar, series')

    rdonwestID = netcdf.defVar(nc_bndry,'RDON_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdonwestID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdonwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdonwestID,'field','RDON, scalar, series')

    rdopwestID = netcdf.defVar(nc_bndry,'RDOP_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdopwestID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rdopwestID,'field','RDOP, scalar, series')

    TAwestID = netcdf.defVar(nc_bndry,'TA_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,TAwestID,'long_name','TOTAL ALKALINITY');
    netcdf.putAtt(nc_bndry,TAwestID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,TAwestID,'field','TA, scalar, series')

    DICwestID = netcdf.defVar(nc_bndry,'DIC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,DICwestID,'long_name','DISSOLVED INORGANIC CARBON');
    netcdf.putAtt(nc_bndry,DICwestID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,DICwestID,'field','DIC, scalar, series')

    CACO3westID = netcdf.defVar(nc_bndry,'CACO3_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,CACO3westID,'long_name','CALCIUM CARBONATE');
    netcdf.putAtt(nc_bndry,CACO3westID,'units','mmol kg-1');
    netcdf.putAtt(nc_bndry,CACO3westID,'field','CACO3, scalar, series')
    
    rpocwestID = netcdf.defVar(nc_bndry,'RPOC_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpocwestID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpocwestID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpocwestID,'field','RPOC, scalar, series')
    
    rponwestID = netcdf.defVar(nc_bndry,'RPON_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rponwestID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponwestID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponwestID,'field','RPON, scalar, series')
    
    rpopwestID = netcdf.defVar(nc_bndry,'RPOP_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpopwestID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopwestID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopwestID,'field','RPOP, scalar, series')
    
    salwestID = netcdf.defVar(nc_bndry,'SAL_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,salwestID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,salwestID,'units','psu');
    netcdf.putAtt(nc_bndry,salwestID,'field','SAL, scalar, series')
    
    sitwestID = netcdf.defVar(nc_bndry,'SIT_west','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,sitwestID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,sitwestID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,sitwestID,'field','SIT, scalar, series')
end
%EAST
if(dir_flag(4)==1)
    aceastID = netcdf.defVar(nc_bndry,'A_C_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,aceastID,'long_name','MX model prey');
    netcdf.putAtt(nc_bndry,aceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,aceastID,'field','A_C, scalar, series');

    achlceastID = netcdf.defVar(nc_bndry,'A_CHLC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,achlceastID,'long_name','MX model CHLC/C of Prey');
    netcdf.putAtt(nc_bndry,achlceastID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,achlceastID,'field','A_CHLC, scalar, series');

    anceastID = netcdf.defVar(nc_bndry,'A_NC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,anceastID,'long_name','MX model N:C of Prey');
    netcdf.putAtt(nc_bndry,anceastID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,anceastID,'field','A_NC, scalar, series');

    apceastID = netcdf.defVar(nc_bndry,'A_PC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,apceastID,'long_name','MX model P:C of Prey');
    netcdf.putAtt(nc_bndry,apceastID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,apceastID,'field','A_PC, scalar, series');

    bsieastID = netcdf.defVar(nc_bndry,'BSI_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,bsieastID,'long_name','BIOGENIC SILICA');
    netcdf.putAtt(nc_bndry,bsieastID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,bsieastID,'field','BSI, scalar, series');

    doeastID = netcdf.defVar(nc_bndry,'DO_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,doeastID,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(nc_bndry,doeastID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,doeastID,'field','DO, scalar, series');

    exdoceastID = netcdf.defVar(nc_bndry,'EXDOC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,exdoceastID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
    netcdf.putAtt(nc_bndry,exdoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,exdoceastID,'field','EXDOC, scalar, series');

    ldoceastID = netcdf.defVar(nc_bndry,'LDOC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldoceastID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
    netcdf.putAtt(nc_bndry,ldoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,ldoceastID,'field','LDOC, scalar, series');

    ldoneastID = netcdf.defVar(nc_bndry,'LDON_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldoneastID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
    netcdf.putAtt(nc_bndry,ldoneastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,ldoneastID,'field','LDON, scalar, series');

    ldopeastID = netcdf.defVar(nc_bndry,'LDOP_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,ldopeastID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
    netcdf.putAtt(nc_bndry,ldopeastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,ldopeastID,'field','LDOP, scalar, series');

    lpoceastID = netcdf.defVar(nc_bndry,'LPOC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpoceastID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
    netcdf.putAtt(nc_bndry,lpoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,lpoceastID,'field','LPOC, scalar, series');

    lponeastID = netcdf.defVar(nc_bndry,'LPON_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lponeastID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
    netcdf.putAtt(nc_bndry,lponeastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,lponeastID,'field','LPON, scalar, series');

    lpopeastID = netcdf.defVar(nc_bndry,'LPOP_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,lpopeastID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
    netcdf.putAtt(nc_bndry,lpopeastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,lpopeastID,'field','LPOP, scalar, series');

    mavgueastID = netcdf.defVar(nc_bndry,'M_AVGCU_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mavgueastID,'long_name','MX model mixotroph average growth rate');
    netcdf.putAtt(nc_bndry,mavgueastID,'units','C/C/d');
    netcdf.putAtt(nc_bndry,mavgueastID,'field','MAVGCU, scalar, series');

    mceastID = netcdf.defVar(nc_bndry,'M_C_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mceastID,'long_name','MX model core mixotroph body biomass carbon');
    netcdf.putAtt(nc_bndry,mceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,mceastID,'field','MC, scalar, series');

    mchlceastID = netcdf.defVar(nc_bndry,'M_CHLC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mchlceastID,'long_name','MX model mixotroph core Chla:C');
    netcdf.putAtt(nc_bndry,mchlceastID,'units','gChl/gC');
    netcdf.putAtt(nc_bndry,mchlceastID,'field','MCHLC, scalar, series');

    mfceastID = netcdf.defVar(nc_bndry,'M_FC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfceastID,'long_name','MX model mixotroph gut relative to body biomass');
    netcdf.putAtt(nc_bndry,mfceastID,'units','gC/gC');
    netcdf.putAtt(nc_bndry,mfceastID,'field','MFC, scalar, series');

    mfchlceastID = netcdf.defVar(nc_bndry,'M_FCHLC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfchlceastID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
    netcdf.putAtt(nc_bndry,mfchlceastID,'units','gChl in the gut /gC');
    netcdf.putAtt(nc_bndry,mfchlceastID,'field','MFCHLC, scalar, series');

    mnceastID = netcdf.defVar(nc_bndry,'M_NC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mnceastID,'long_name','MX model mixotroph N:C');
    netcdf.putAtt(nc_bndry,mnceastID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mnceastID,'field','MNC, scalar, series');

    mpceastID = netcdf.defVar(nc_bndry,'M_PC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mpceastID,'long_name','MX model mixotroph P:C');
    netcdf.putAtt(nc_bndry,mpceastID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mpceastID,'field','MPC, scalar, series');

    mfnceastID = netcdf.defVar(nc_bndry,'MFNC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfnceastID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfnceastID,'units','gN/gC');
    netcdf.putAtt(nc_bndry,mfnceastID,'field','MFNC, scalar, series');

    mfpceastID = netcdf.defVar(nc_bndry,'MFPC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,mfpceastID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
    netcdf.putAtt(nc_bndry,mfpceastID,'units','gP/gC');
    netcdf.putAtt(nc_bndry,mfpceastID,'field','MFPC, scalar, series');

    nh4teastID = netcdf.defVar(nc_bndry,'NH4T_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,nh4teastID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
    netcdf.putAtt(nc_bndry,nh4teastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,nh4teastID,'field','NH4T, scalar, series');

    no23eastID = netcdf.defVar(nc_bndry,'NO23_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,no23eastID,'long_name','NITRITE + NITRATE');
    netcdf.putAtt(nc_bndry,no23eastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,no23eastID,'field','NO23, scalar, series')

    o2eqeastID = netcdf.defVar(nc_bndry,'O2EQ_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,o2eqeastID,'long_name','AQUEOUS SOD');
    netcdf.putAtt(nc_bndry,o2eqeastID,'units','mg O2 L-1');
    netcdf.putAtt(nc_bndry,o2eqeastID,'field','O2EQ, scalar, series')

    phyt1eastID = netcdf.defVar(nc_bndry,'PHYT1_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt1eastID,'long_name','WINTER DIATOMS (PHYT1)');
    netcdf.putAtt(nc_bndry,phyt1eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt1eastID,'field','PHYT1, scalar, series')

    phyt2eastID = netcdf.defVar(nc_bndry,'PHYT2_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt2eastID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
    netcdf.putAtt(nc_bndry,phyt2eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt2eastID,'field','PHYT2, scalar, series')

    phyt3eastID = netcdf.defVar(nc_bndry,'PHYT3_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,phyt3eastID,'long_name','FALL ASSEMBLAGE (PHYT3)');
    netcdf.putAtt(nc_bndry,phyt3eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,phyt3eastID,'field','PHYT3, scalar, series')
 
    zoo1eastID = netcdf.defVar(nc_bndry,'ZOO1_east','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo1eastID,'long_name','SMALL ZOOPLANKTON (ZOO1)');
    netcdf.putAtt(nc_bndry,zoo1eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo1eastID,'field','ZOO1, scalar, series')

    zoo2eastID = netcdf.defVar(nc_bndry,'ZOO2_east','float',[xrhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,zoo2eastID,'long_name','LARGE ZOOPLANKTON (ZOO2)');
    netcdf.putAtt(nc_bndry,zoo2eastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,zoo2eastID,'field','ZOO2, scalar, series')

    po4teastID = netcdf.defVar(nc_bndry,'PO4T_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,po4teastID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
    netcdf.putAtt(nc_bndry,po4teastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,po4teastID,'field','PO4T, scalar, series')

    rdoceastID = netcdf.defVar(nc_bndry,'RDOC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdoceastID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
    netcdf.putAtt(nc_bndry,rdoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rdoceastID,'field','RDOC, scalar, series')

    rdoneastID = netcdf.defVar(nc_bndry,'RDON_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdoneastID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
    netcdf.putAtt(nc_bndry,rdoneastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdoneastID,'field','RDON, scalar, series')

    rdopeastID = netcdf.defVar(nc_bndry,'RDOP_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rdopeastID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
    netcdf.putAtt(nc_bndry,rdopeastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rdopeastID,'field','RDOP, scalar, series')

    TAeastID = netcdf.defVar(nc_bndry,'TA_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,TAeastID,'long_name','TOTAL ALKALINITY');
    netcdf.putAtt(nc_bndry,TAeastID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,TAeastID,'field','TA, scalar, series')

    DICeastID = netcdf.defVar(nc_bndry,'DIC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,DICeastID,'long_name','DISSOLVED INORGANIC CARBON');
    netcdf.putAtt(nc_bndry,DICeastID,'units','umol L-1');
    netcdf.putAtt(nc_bndry,DICeastID,'field','DIC, scalar, series')

    CACO3eastID = netcdf.defVar(nc_bndry,'CACO3_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,CACO3eastID,'long_name','CALCIUM CARBONATE');
    netcdf.putAtt(nc_bndry,CACO3eastID,'units','mmol kg-1');
    netcdf.putAtt(nc_bndry,CACO3eastID,'field','CACO3, scalar, series')
    
    rpoceastID = netcdf.defVar(nc_bndry,'RPOC_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpoceastID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
    netcdf.putAtt(nc_bndry,rpoceastID,'units','mg C L-1');
    netcdf.putAtt(nc_bndry,rpoceastID,'field','RPOC, scalar, series')
    
    rponeastID = netcdf.defVar(nc_bndry,'RPON_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rponeastID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
    netcdf.putAtt(nc_bndry,rponeastID,'units','mg N L-1');
    netcdf.putAtt(nc_bndry,rponeastID,'field','RPON, scalar, series')
    
    rpopeastID = netcdf.defVar(nc_bndry,'RPOP_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,rpopeastID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
    netcdf.putAtt(nc_bndry,rpopeastID,'units','mg P L-1');
    netcdf.putAtt(nc_bndry,rpopeastID,'field','RPOP, scalar, series')
    
    saleastID = netcdf.defVar(nc_bndry,'SAL_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,saleastID,'long_name','SALINITY');
    netcdf.putAtt(nc_bndry,saleastID,'units','psu');
    netcdf.putAtt(nc_bndry,saleastID,'field','SAL, scalar, series')
    
    siteastID = netcdf.defVar(nc_bndry,'SIT_east','float',[erhodimID s_rhodimID biotdimID]);
    netcdf.putAtt(nc_bndry,siteastID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
    netcdf.putAtt(nc_bndry,siteastID,'units','mg SI L-1');
    netcdf.putAtt(nc_bndry,siteastID,'field','SIT, scalar, series')
end

netcdf.close(nc_bndry)

end

