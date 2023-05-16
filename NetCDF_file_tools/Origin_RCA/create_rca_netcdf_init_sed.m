function create_rca_netcdf_init_sed(fn,gn,t)
% Generate a blank netcdf file containing RCA Sediment initial conditions.
% Yuren Chen, 2021/06/22, Guangzhou
%
% Get some grid info. 
  [LP,MP]=size(gn.lon_rho); 
  NX  = LP+2;
  NY = MP+2; 
  NZ    = gn.N;
  NZZ = NZ+1;
  group = 3;
  NSYS = 26;
  
%% create bndry file
nc_init=netcdf.create(fn,'clobber');
if isempty(nc_init), return, end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'history', ['Created by cyr on ' datestr(now)]);
netcdf.putAtt(nc_init,netcdf.getConstant('NC_GLOBAL'),'type', 'Nutrients data from Florida Water Atlas');
%% Dimensions:

disp(' ## Defining Dimensions...')
 
nxdimID = netcdf.defDim(nc_init,'NX',NX);
nydimID = netcdf.defDim(nc_init,'NY',NY);
nzdimID = netcdf.defDim(nc_init,'NZ',NZ);
nzzdimID = netcdf.defDim(nc_init,'NZZ',NZZ);
groupdimID = netcdf.defDim(nc_init,'group',group);
nsysdimID = netcdf.defDim(nc_init,'NSYS',NSYS);
timedimID = netcdf.defDim(nc_init,'TIME',t);

%% Variables and attributes:

TIMEID = netcdf.defVar(nc_init,'TIME','float',[timedimID]);
netcdf.putAtt(nc_init,TIMEID,'long_name','time since initialization');
netcdf.putAtt(nc_init,TIMEID,'units','days');
netcdf.putAtt(nc_init,TIMEID,'field','TIME, scalar, series');

BFORMAXSID = netcdf.defVar(nc_init,'BFORMAXS','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,BFORMAXSID,'long_name',' ');
netcdf.putAtt(nc_init,BFORMAXSID,'units',' ');
netcdf.putAtt(nc_init,BFORMAXSID,'field','BFORMAXS, scalar, series');

BNTHSTR1SID = netcdf.defVar(nc_init,'BNTHSTR1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,BNTHSTR1SID,'long_name','benthic stress');
netcdf.putAtt(nc_init,BNTHSTR1SID,'units',' ');
netcdf.putAtt(nc_init,BNTHSTR1SID,'field','BNTHSTR1S, scalar, series');

CH41TM1SID = netcdf.defVar(nc_init,'CH41TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,CH41TM1SID,'long_name','methane');
netcdf.putAtt(nc_init,CH41TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,CH41TM1SID,'field','CH41TM1S, scalar, series');

CH4T1TM1SID = netcdf.defVar(nc_init,'CH4T1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,CH4T1TM1SID,'long_name','dissolved CH4');
netcdf.putAtt(nc_init,CH4T1TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,CH4T1TM1SID,'field','CH4T1TM1S, scalar, series');

CH4T2TM1SID = netcdf.defVar(nc_init,'CH4T2TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,CH4T2TM1SID,'long_name','dissolved CH4');
netcdf.putAtt(nc_init,CH4T2TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,CH4T2TM1SID,'field','CH4T2TM1S, scalar, series');

CPOCID = netcdf.defVar(nc_init,'CPOC','float',[nxdimID nydimID groupdimID timedimID]);
netcdf.putAtt(nc_init,CPOCID,'long_name','particulate organic C');
netcdf.putAtt(nc_init,CPOCID,'units','MG C/M^3');
netcdf.putAtt(nc_init,CPOCID,'field','CPOC, scalar, series');

CPONID = netcdf.defVar(nc_init,'CPON','float',[nxdimID nydimID groupdimID timedimID]);
netcdf.putAtt(nc_init,CPONID,'long_name','particulate organic N');
netcdf.putAtt(nc_init,CPONID,'units','MG N/M^3');
netcdf.putAtt(nc_init,CPONID,'field','CPON, scalar, series');

CPOPID = netcdf.defVar(nc_init,'CPOP','float',[nxdimID nydimID groupdimID timedimID]);
netcdf.putAtt(nc_init,CPOPID,'long_name','particulate organic P');
netcdf.putAtt(nc_init,CPOPID,'units','MG P/M^3');
netcdf.putAtt(nc_init,CPOPID,'field','CPOP, scalar, series');

CPOSID = netcdf.defVar(nc_init,'CPOS','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,CPOSID,'long_name','biogenic Si');
netcdf.putAtt(nc_init,CPOSID,'units','MG Si/M^3');
netcdf.putAtt(nc_init,CPOSID,'field','CPOS, scalar, series');

CTEMPID = netcdf.defVar(nc_init,'CTEMP','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,CTEMPID,'long_name','temperature');
netcdf.putAtt(nc_init,CTEMPID,'units','Deg C');
netcdf.putAtt(nc_init,CTEMPID,'field','CTEMP, scalar, series');

DD0TM1SID = netcdf.defVar(nc_init,'DD0TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,DD0TM1SID,'long_name','');
netcdf.putAtt(nc_init,DD0TM1SID,'units','');
netcdf.putAtt(nc_init,DD0TM1SID,'field','DD0TM1S, scalar, series');

HS1TM1SID = netcdf.defVar(nc_init,'HS1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,HS1TM1SID,'long_name','dissolved H2S');
netcdf.putAtt(nc_init,HS1TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,HS1TM1SID,'field','HS1TM1S, scalar, series');

HST1TM1SID = netcdf.defVar(nc_init,'HST1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,HST1TM1SID,'long_name','dissolved H2S');
netcdf.putAtt(nc_init,HST1TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,HST1TM1SID,'field','HST1TM1S, scalar, series');

HST2TM1SID = netcdf.defVar(nc_init,'HST2TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,HST2TM1SID,'long_name','dissolved H2S');
netcdf.putAtt(nc_init,HST2TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,HST2TM1SID,'field','HST2TM1S, scalar, series');

ISWBNTHSID = netcdf.defVar(nc_init,'ISWBNTHS','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,ISWBNTHSID,'long_name','');
netcdf.putAtt(nc_init,ISWBNTHSID,'units','');
netcdf.putAtt(nc_init,ISWBNTHSID,'field','ISWBNTHS, scalar, series')

NH41TM1SID = netcdf.defVar(nc_init,'NH41TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,NH41TM1SID,'long_name','dissolved NH4');
netcdf.putAtt(nc_init,NH41TM1SID,'units','MG N/M^3');
netcdf.putAtt(nc_init,NH41TM1SID,'field','NH41TM1S, scalar, series')

NH4T1TM1SID = netcdf.defVar(nc_init,'NH4T1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,NH4T1TM1SID,'long_name','dissolved NH4');
netcdf.putAtt(nc_init,NH4T1TM1SID,'units','MG N/M^3');
netcdf.putAtt(nc_init,NH4T1TM1SID,'field','NH4T1TM1S, scalar, series')

NH4T2TM1SID = netcdf.defVar(nc_init,'NH4T2TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,NH4T2TM1SID,'long_name','dissolved NH4');
netcdf.putAtt(nc_init,NH4T2TM1SID,'units','MG N/M^3');
netcdf.putAtt(nc_init,NH4T2TM1SID,'field','NH4T2TM1S, scalar, series')

NO31TM1SID = netcdf.defVar(nc_init,'NO31TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,NO31TM1SID,'long_name','dissolved NO3');
netcdf.putAtt(nc_init,NO31TM1SID,'units','MG N/M^3');
netcdf.putAtt(nc_init,NO31TM1SID,'field','NO31TM1S, scalar, series')

NO3T1TM1SID = netcdf.defVar(nc_init,'NO3T1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,NO3T1TM1SID,'long_name','dissolved NO3');
netcdf.putAtt(nc_init,NO3T1TM1SID,'units','MG N/M^3');
netcdf.putAtt(nc_init,NO3T1TM1SID,'field','NO3T1TM1S, scalar, series')

NO3T2TM1SID = netcdf.defVar(nc_init,'NO3T2TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,NO3T2TM1SID,'long_name','dissolved NO3');
netcdf.putAtt(nc_init,NO3T2TM1SID,'units','MG N/M^3');
netcdf.putAtt(nc_init,NO3T2TM1SID,'field','NO3T2TM1S, scalar, series')

O20TM1SID = netcdf.defVar(nc_init,'O20TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,O20TM1SID,'long_name','dissolved O2');
netcdf.putAtt(nc_init,O20TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,O20TM1SID,'field','O20TM1S, scalar, series')

PO41TM1SID = netcdf.defVar(nc_init,'PO41TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,PO41TM1SID,'long_name','dissolved PO4');
netcdf.putAtt(nc_init,PO41TM1SID,'units','MG P/M^3');
netcdf.putAtt(nc_init,PO41TM1SID,'field','PO41TM1S, scalar, series')

PO4T1TM1SID = netcdf.defVar(nc_init,'PO4T1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,PO4T1TM1SID,'long_name','dissolved PO4');
netcdf.putAtt(nc_init,PO4T1TM1SID,'units','MG P/M^3');
netcdf.putAtt(nc_init,PO4T1TM1SID,'field','PO4T1TM1S, scalar, series')

PO4T2TM1SID = netcdf.defVar(nc_init,'PO4T2TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,PO4T2TM1SID,'long_name','dissolved PO4');
netcdf.putAtt(nc_init,PO4T2TM1SID,'units','MG P/M^3');
netcdf.putAtt(nc_init,PO4T2TM1SID,'field','PO4T2TM1S, scalar, series')

SI1TM1SID = netcdf.defVar(nc_init,'SI1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,SI1TM1SID,'long_name','dissolved SI');
netcdf.putAtt(nc_init,SI1TM1SID,'units','MG SI/M^3');
netcdf.putAtt(nc_init,SI1TM1SID,'field','SI1TM1S, scalar, series')

SIT1TM1SID = netcdf.defVar(nc_init,'SIT1TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,SIT1TM1SID,'long_name','dissolved SI');
netcdf.putAtt(nc_init,SIT1TM1SID,'units','MG SI/M^3');
netcdf.putAtt(nc_init,SIT1TM1SID,'field','SIT1TM1S, scalar, series')

SIT2TM1SID = netcdf.defVar(nc_init,'SIT2TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,SIT2TM1SID,'long_name','dissolved SI');
netcdf.putAtt(nc_init,SIT2TM1SID,'units','MG SI/M^3');
netcdf.putAtt(nc_init,SIT2TM1SID,'field','SIT2TM1S, scalar, series')

SO4T2TM1SID = netcdf.defVar(nc_init,'SO4T2TM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,SO4T2TM1SID,'long_name','sulfate');
netcdf.putAtt(nc_init,SO4T2TM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,SO4T2TM1SID,'field','SO4T2TM1S, scalar, series')

SODTM1SID = netcdf.defVar(nc_init,'SODTM1S','float',[nxdimID nydimID timedimID]);
netcdf.putAtt(nc_init,SODTM1SID,'long_name','sediment oxygen demand');
netcdf.putAtt(nc_init,SODTM1SID,'units','MG O2/M^3');
netcdf.putAtt(nc_init,SODTM1SID,'field','SODTM1S, scalar, series')

netcdf.close(nc_init)
