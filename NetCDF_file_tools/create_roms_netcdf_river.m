function create_roms_netcdf_river(river_file,gn,num_rivers,num_river_times)

%get some grid info
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
N  = gn.N;

%create init file
nc_river=netcdf.create(river_file, 'clobber');
if isempty(nc_river)
    disp([' ## Unable to create ROMS Rivers NetCDF file.'])
    return
end

%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_river,netcdf.getConstant('NC_GLOBAL'),'type','ROMS Rivers forcing file');
netcdf.putAtt(nc_river,netcdf.getConstant('NC_GLOBAL'),'history',['Created by ', mfilename ', on ', datestr(now)]);
netcdf.putAtt(nc_river,netcdf.getConstant('NC_GLOBAL'),'title','ROMS Application')

%% Dimensions:
disp(' ## Defining Dimensions...')

xpsidimID = netcdf.defDim(nc_river,'xpsi',L);
xrhodimID = netcdf.defDim(nc_river,'xrho',LP);
xudimID   = netcdf.defDim(nc_river,'xu',L);
xvdimID   = netcdf.defDim(nc_river,'xv',LP);

epsidimID = netcdf.defDim(nc_river,'epsi',M);
erhodimID = netcdf.defDim(nc_river,'erho',MP);
eudimID   = netcdf.defDim(nc_river,'eu',MP);
evdimID   = netcdf.defDim(nc_river,'ev',M);

s_rhodimID = netcdf.defDim(nc_river,'s_rho',N);
s_wdimID = netcdf.defDim(nc_river,'s_w',N+1);
riverdimID = netcdf.defDim(nc_river,'river',num_rivers);
river_timedimID = netcdf.defDim(nc_river,'river_time',num_river_times);
onedimID = netcdf.defDim(nc_river,'one',1);

%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')

theta_sID = netcdf.defVar(nc_river,'theta_s','double',onedimID);
netcdf.putAtt(nc_river,theta_sID,'long_name','S-coordinate surface control parameter');
netcdf.putAtt(nc_river,theta_sID,'units','1');

theta_bID = netcdf.defVar(nc_river,'theta_b','double',onedimID);
netcdf.putAtt(nc_river,theta_bID,'long_name','S-coordinate bottom control parameter');
netcdf.putAtt(nc_river,theta_bID,'units','1');

tcline_ID = netcdf.defVar(nc_river,'Tcline','double',onedimID);
netcdf.putAtt(nc_river,tcline_ID,'long_name','S-coordinate surface/bottom layer width');
netcdf.putAtt(nc_river,tcline_ID,'units','meter');

hc_ID = netcdf.defVar(nc_river,'hc','double',onedimID);
netcdf.putAtt(nc_river,hc_ID,'long_name','S-coordinate parameter, critical depth');
netcdf.putAtt(nc_river,hc_ID,'units','meter');

Cs_rID = netcdf.defVar(nc_river,'Cs_r','double',s_rhodimID);
netcdf.putAtt(nc_river,Cs_rID,'long_name','S-coordinate stretching curves at RHO-points');
netcdf.putAtt(nc_river,Cs_rID,'units','1');
netcdf.putAtt(nc_river,Cs_rID,'valid_min',-1);
netcdf.putAtt(nc_river,Cs_rID,'valid_max',0);
netcdf.putAtt(nc_river,Cs_rID,'field','Cs_r, scalar');

Cs_wID = netcdf.defVar(nc_river,'Cs_w','double',s_wdimID);
netcdf.putAtt(nc_river,Cs_wID,'long_name','S-coordinate stretching curves at W-points');
netcdf.putAtt(nc_river,Cs_wID,'units','1');
netcdf.putAtt(nc_river,Cs_wID,'valid_min',-1);
netcdf.putAtt(nc_river,Cs_wID,'valid_max',0);
netcdf.putAtt(nc_river,Cs_wID,'field','Cs_w, scalar');

sc_rID = netcdf.defVar(nc_river,'sc_r','double',s_rhodimID);
netcdf.putAtt(nc_river,sc_rID,'long_name','S-coordinate at RHO-points');
netcdf.putAtt(nc_river,sc_rID,'units','1');
netcdf.putAtt(nc_river,sc_rID,'valid_min',-1);
netcdf.putAtt(nc_river,sc_rID,'valid_max',0);
netcdf.putAtt(nc_river,sc_rID,'field','sc_r, scalar');

sc_wID = netcdf.defVar(nc_river,'sc_w','double',s_wdimID);
netcdf.putAtt(nc_river,sc_wID,'long_name','S-coordinate at W-points');
netcdf.putAtt(nc_river,sc_wID,'units','1');
netcdf.putAtt(nc_river,sc_wID,'valid_min',-1);
netcdf.putAtt(nc_river,sc_wID,'valid_max',0);
netcdf.putAtt(nc_river,sc_wID,'field','sc_w, scalar');

river_ID = netcdf.defVar(nc_river,'river','double',riverdimID);
netcdf.putAtt(nc_river,river_ID,'long_name','river_runoff identification number');
netcdf.putAtt(nc_river,river_ID,'units','nondimensional');
netcdf.putAtt(nc_river,river_ID,'field','num_rivers, scalar, series');

river_timeID = netcdf.defVar(nc_river,'river_time','double',river_timedimID);
netcdf.putAtt(nc_river,river_timeID,'long_name','river time');
netcdf.putAtt(nc_river,river_timeID,'units','days');
netcdf.putAtt(nc_river,river_timeID,'field','river_time, scalar, series');

river_XpositionID = netcdf.defVar(nc_river,'river_Xposition','double',riverdimID);
netcdf.putAtt(nc_river,river_XpositionID,'long_name','river runoff  XI-positions at RHO-points');
netcdf.putAtt(nc_river,river_XpositionID,'units','scalar');
netcdf.putAtt(nc_river,river_XpositionID,'time','river_time');
netcdf.putAtt(nc_river,river_XpositionID,'field','river runoff XI position, scalar, series');

river_EpositionID = netcdf.defVar(nc_river,'river_Eposition','double',riverdimID);
netcdf.putAtt(nc_river,river_EpositionID,'long_name','river runoff  ETA-positions at RHO-points');
netcdf.putAtt(nc_river,river_EpositionID,'units','scalar');
netcdf.putAtt(nc_river,river_EpositionID,'time','river_time');
netcdf.putAtt(nc_river,river_EpositionID,'field','river runoff ETA position, scalar, series');

river_directionID = netcdf.defVar(nc_river,'river_direction','double',riverdimID);
netcdf.putAtt(nc_river,river_directionID,'long_name','river runoff direction, XI=0, ETA>0');
netcdf.putAtt(nc_river,river_directionID,'units','scalar');
netcdf.putAtt(nc_river,river_directionID,'time','river_time');
netcdf.putAtt(nc_river,river_directionID,'field','river runoff direction, scalar, series');

river_VshapeID = netcdf.defVar(nc_river,'river_Vshape','double',[riverdimID s_rhodimID]);
netcdf.putAtt(nc_river,river_VshapeID,'long_name','river runoff mass transport vertical profile');
netcdf.putAtt(nc_river,river_VshapeID,'units','scalar');
netcdf.putAtt(nc_river,river_VshapeID,'time','river_time');
netcdf.putAtt(nc_river,river_VshapeID,'field','river runoff vertical profile, scalar, series');

river_transportID = netcdf.defVar(nc_river,'river_transport','double',[riverdimID river_timedimID]);
netcdf.putAtt(nc_river,river_transportID,'long_name','river runoff mass transport');
netcdf.putAtt(nc_river,river_transportID,'units','meter^3/s');
netcdf.putAtt(nc_river,river_transportID,'time','river_time');
netcdf.putAtt(nc_river,river_transportID,'field','river runoff mass transport, scalar, series');

acriverID = netcdf.defVar(nc_river,'river_A_C','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,acriverID,'long_name','MX model prey');
netcdf.putAtt(nc_river,acriverID,'units','mg C L-1');
netcdf.putAtt(nc_river, acriverID,'time','river_time');
netcdf.putAtt(nc_river,acriverID,'field','river_A_C, scalar, series');

achlcriverID = netcdf.defVar(nc_river,'river_A_CHLC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,achlcriverID,'long_name','MX model CHLC/C of Prey');
netcdf.putAtt(nc_river,achlcriverID,'units','gChl/gC');
netcdf.putAtt(nc_river, achlcriverID,'time','river_time');
netcdf.putAtt(nc_river,achlcriverID,'field','river_A_CHLC, scalar, series');

ancriverID = netcdf.defVar(nc_river,'river_A_NC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,ancriverID,'long_name','MX model N:C of Prey');
netcdf.putAtt(nc_river,ancriverID,'units','gN/gC');
netcdf.putAtt(nc_river, ancriverID,'time','river_time');
netcdf.putAtt(nc_river,ancriverID,'field','river_A_NC, scalar, series');

apcriverID = netcdf.defVar(nc_river,'river_A_PC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,apcriverID,'long_name','MX model P:C of Prey');
netcdf.putAtt(nc_river,apcriverID,'units','gP/gC');
netcdf.putAtt(nc_river, apcriverID,'time','river_time');
netcdf.putAtt(nc_river,apcriverID,'field','river_A_PC, scalar, series');

bsiriverID = netcdf.defVar(nc_river,'river_BSI','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,bsiriverID,'long_name','BIOGENIC SILICA');
netcdf.putAtt(nc_river,bsiriverID,'units','mg SI L-1');
netcdf.putAtt(nc_river, bsiriverID,'time','river_time');
netcdf.putAtt(nc_river,bsiriverID,'field','river_BSI, scalar, series');

doriverID = netcdf.defVar(nc_river,'river_DO','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,doriverID,'long_name','DISSOLVED OXYGEN');
netcdf.putAtt(nc_river,doriverID,'units','mg O2 L-1');
netcdf.putAtt(nc_river, doriverID,'time','river_time');
netcdf.putAtt(nc_river,doriverID,'field','river_DO, scalar, series');

exdocriverID = netcdf.defVar(nc_river,'river_EXDOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,exdocriverID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
netcdf.putAtt(nc_river,exdocriverID,'units','mg C L-1');
netcdf.putAtt(nc_river, exdocriverID,'time','river_time');
netcdf.putAtt(nc_river,exdocriverID,'field','river_exdoc, scalar, series');

ldocriverID = netcdf.defVar(nc_river,'river_LDOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,ldocriverID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
netcdf.putAtt(nc_river,ldocriverID,'units','mg C L-1');
netcdf.putAtt(nc_river, ldocriverID,'time','river_time');
netcdf.putAtt(nc_river,ldocriverID,'field','river_ldoc, scalar, series');

ldonriverID = netcdf.defVar(nc_river,'river_LDON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,ldonriverID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
netcdf.putAtt(nc_river,ldonriverID,'units','mg N L-1');
netcdf.putAtt(nc_river, ldonriverID,'time','river_time');
netcdf.putAtt(nc_river,ldonriverID,'field','river_ldon, scalar, series');

ldopriverID = netcdf.defVar(nc_river,'river_LDOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,ldopriverID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
netcdf.putAtt(nc_river,ldopriverID,'units','mg P L-1');
netcdf.putAtt(nc_river, ldopriverID,'time','river_time');
netcdf.putAtt(nc_river,ldopriverID,'field','river_ldop, scalar, series');

lpocriverID = netcdf.defVar(nc_river,'river_LPOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,lpocriverID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
netcdf.putAtt(nc_river,lpocriverID,'units','mg C L-1');
netcdf.putAtt(nc_river, lpocriverID,'time','river_time');
netcdf.putAtt(nc_river,lpocriverID,'field','river_lpoc, scalar, series');

lponriverID = netcdf.defVar(nc_river,'river_LPON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,lponriverID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
netcdf.putAtt(nc_river,lponriverID,'units','mg N L-1');
netcdf.putAtt(nc_river, lponriverID,'time','river_time');
netcdf.putAtt(nc_river,lponriverID,'field','river_lpon, scalar, series');

lpopriverID = netcdf.defVar(nc_river,'river_LPOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,lpopriverID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
netcdf.putAtt(nc_river,lpopriverID,'units','mg P L-1');
netcdf.putAtt(nc_river, lpopriverID,'time','river_time');
netcdf.putAtt(nc_river,lpopriverID,'field','river_lpop, scalar, series');

mavguriverID = netcdf.defVar(nc_river,'river_M_AVGCU','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mavguriverID,'long_name','MX model mixotroph average growth rate');
netcdf.putAtt(nc_river,mavguriverID,'units','C/C/d');
netcdf.putAtt(nc_river, mavguriverID,'time','river_time');
netcdf.putAtt(nc_river,mavguriverID,'field','river_mavgu, scalar, series');

mcriverID = netcdf.defVar(nc_river,'river_M_C','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mcriverID,'long_name','MX model core mixotroph body biomass carbon');
netcdf.putAtt(nc_river,mcriverID,'units','mg C L-1');
netcdf.putAtt(nc_river, mcriverID,'time','river_time');
netcdf.putAtt(nc_river,mcriverID,'field','river_mc, scalar, series');

mchlcriverID = netcdf.defVar(nc_river,'river_M_CHLC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mchlcriverID,'long_name','MX model mixotroph core Chla:C');
netcdf.putAtt(nc_river,mchlcriverID,'units','gChl/gC');
netcdf.putAtt(nc_river, mchlcriverID,'time','river_time');
netcdf.putAtt(nc_river,mchlcriverID,'field','river_mchlc, scalar, series');

mfcriverID = netcdf.defVar(nc_river,'river_M_FC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mfcriverID,'long_name','MX model mixotroph gut relative to body biomass');
netcdf.putAtt(nc_river,mfcriverID,'units','gC/gC');
netcdf.putAtt(nc_river, mfcriverID,'time','river_time');
netcdf.putAtt(nc_river,mfcriverID,'field','river_mfc, scalar, series');

mfchlcriverID = netcdf.defVar(nc_river,'river_M_FCHLC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mfchlcriverID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
netcdf.putAtt(nc_river,mfchlcriverID,'units','gChl in the gut /gC');
netcdf.putAtt(nc_river, mfchlcriverID,'time','river_time');
netcdf.putAtt(nc_river,mfchlcriverID,'field','river_mfchlc, scalar, series');

mncriverID = netcdf.defVar(nc_river,'river_M_NC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mncriverID,'long_name','MX model mixotroph N:C');
netcdf.putAtt(nc_river,mncriverID,'units','gN/gC');
netcdf.putAtt(nc_river, mncriverID,'time','river_time');
netcdf.putAtt(nc_river,mncriverID,'field','river_mnc, scalar, series');

mpcriverID = netcdf.defVar(nc_river,'river_M_PC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mpcriverID,'long_name','MX model mixotroph P:C');
netcdf.putAtt(nc_river,mpcriverID,'units','gP/gC');
netcdf.putAtt(nc_river, mpcriverID,'time','river_time');
netcdf.putAtt(nc_river,mpcriverID,'field','river_mpc, scalar, series');

mriver_filecriverID = netcdf.defVar(nc_river,'river_MFNC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mriver_filecriverID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
netcdf.putAtt(nc_river,mriver_filecriverID,'units','gN/gC');
netcdf.putAtt(nc_river, mriver_filecriverID,'time','river_time');
netcdf.putAtt(nc_river,mriver_filecriverID,'field','river_mriver_filec, scalar, series');

mfpcriverID = netcdf.defVar(nc_river,'river_MFPC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,mfpcriverID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
netcdf.putAtt(nc_river,mfpcriverID,'units','gP/gC');
netcdf.putAtt(nc_river, mfpcriverID,'time','river_time');
netcdf.putAtt(nc_river,mfpcriverID,'field','river_mfpc, scalar, series');

nh4triverID = netcdf.defVar(nc_river,'river_NH4T','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,nh4triverID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
netcdf.putAtt(nc_river,nh4triverID,'units','mg N L-1');
netcdf.putAtt(nc_river, nh4triverID,'time','river_time');
netcdf.putAtt(nc_river,nh4triverID,'field','river_nh4t, scalar, series');

no23riverID = netcdf.defVar(nc_river,'river_NO23','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,no23riverID,'long_name','NITRITE + NITRATE');
netcdf.putAtt(nc_river,no23riverID,'units','mg N L-1');
netcdf.putAtt(nc_river, no23riverID,'time','river_time');
netcdf.putAtt(nc_river,no23riverID,'field','river_no23, scalar, series')

o2eqriverID = netcdf.defVar(nc_river,'river_O2EQ','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,o2eqriverID,'long_name','AQUEOUS SOD');
netcdf.putAtt(nc_river,o2eqriverID,'units','mg O2 L-1');
netcdf.putAtt(nc_river, o2eqriverID,'time','river_time');
netcdf.putAtt(nc_river,o2eqriverID,'field','river_o2eq, scalar, series')

phyt1riverID = netcdf.defVar(nc_river,'river_PHYT1','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,phyt1riverID,'long_name','WINTER DIATOMS (PHYT1)');
netcdf.putAtt(nc_river,phyt1riverID,'units','mg C L-1');
netcdf.putAtt(nc_river, phyt1riverID,'time','river_time');
netcdf.putAtt(nc_river,phyt1riverID,'field','river_phyt1, scalar, series')

phyt2riverID = netcdf.defVar(nc_river,'river_PHYT2','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,phyt2riverID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
netcdf.putAtt(nc_river,phyt2riverID,'units','mg C L-1');
netcdf.putAtt(nc_river, phyt2riverID,'time','river_time');
netcdf.putAtt(nc_river,phyt2riverID,'field','river_phyt2, scalar, series')

phyt3riverID = netcdf.defVar(nc_river,'river_PHYT3','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,phyt3riverID,'long_name','FALL ASSEMBLAGE (PHYT3)');
netcdf.putAtt(nc_river,phyt3riverID,'units','mg C L-1');
netcdf.putAtt(nc_river, phyt3riverID,'time','river_time');
netcdf.putAtt(nc_river,phyt3riverID,'field','river_phyt3, scalar, series')

zoo1riverID = netcdf.defVar(nc_river,'river_ZOO1','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,zoo1riverID,'long_name','SMALL ZOOPLANKTON (ZOO1)');
netcdf.putAtt(nc_river,zoo1riverID,'units','mg C L-1');
netcdf.putAtt(nc_river, zoo1riverID,'time','river_time');
netcdf.putAtt(nc_river,zoo1riverID,'field','river_zoo1, scalar, series')

zoo2riverID = netcdf.defVar(nc_river,'river_ZOO2','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,zoo2riverID,'long_name','LARGE ZOOPLANKTON (ZOO2)');
netcdf.putAtt(nc_river,zoo2riverID,'units','mg C L-1');
netcdf.putAtt(nc_river, zoo2riverID,'time','river_time');
netcdf.putAtt(nc_river,zoo2riverID,'field','river_zoo2, scalar, series')

po4triverID = netcdf.defVar(nc_river,'river_PO4T','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,po4triverID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
netcdf.putAtt(nc_river,po4triverID,'units','mg P L-1');
netcdf.putAtt(nc_river, po4triverID,'time','river_time');
netcdf.putAtt(nc_river,po4triverID,'field','river_po4t, scalar, series')

rdocriverID = netcdf.defVar(nc_river,'river_RDOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,rdocriverID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
netcdf.putAtt(nc_river,rdocriverID,'units','mg C L-1');
netcdf.putAtt(nc_river, rdocriverID,'time','river_time');
netcdf.putAtt(nc_river,rdocriverID,'field','river_rdoc, scalar, series')

rdonriverID = netcdf.defVar(nc_river,'river_RDON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,rdonriverID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
netcdf.putAtt(nc_river,rdonriverID,'units','mg N L-1');
netcdf.putAtt(nc_river, rdonriverID,'time','river_time');
netcdf.putAtt(nc_river,rdonriverID,'field','river_rdon, scalar, series')

rdopriverID = netcdf.defVar(nc_river,'river_RDOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,rdopriverID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
netcdf.putAtt(nc_river,rdopriverID,'units','mg P L-1');
netcdf.putAtt(nc_river, rdopriverID,'time','river_time');
netcdf.putAtt(nc_river,rdopriverID,'field','river_rdop, scalar, series')

TAriverID = netcdf.defVar(nc_river,'river_TA','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,TAriverID,'long_name','TOTAL ALKALINITY');
netcdf.putAtt(nc_river,TAriverID,'units','umol L-1');
netcdf.putAtt(nc_river, TAriverID,'time','river_time');
netcdf.putAtt(nc_river,TAriverID,'field','river_TA, scalar, series')

DICriverID = netcdf.defVar(nc_river,'river_DIC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,DICriverID,'long_name','DISSOLVED INORGANIC CARBON');
netcdf.putAtt(nc_river,DICriverID,'units','umol L-1');
netcdf.putAtt(nc_river, DICriverID,'time','river_time');
netcdf.putAtt(nc_river,DICriverID,'field','river_DIC, scalar, series')

CACO3riverID = netcdf.defVar(nc_river,'river_CACO3','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,CACO3riverID,'long_name','CALCIUM CARBONATE');
netcdf.putAtt(nc_river,CACO3riverID,'units','mmol L-1');
netcdf.putAtt(nc_river, CACO3riverID,'time','river_time');
netcdf.putAtt(nc_river,CACO3riverID,'field','river_CACO3, scalar, series')

rpocriverID = netcdf.defVar(nc_river,'river_RPOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,rpocriverID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
netcdf.putAtt(nc_river,rpocriverID,'units','mg C L-1');
netcdf.putAtt(nc_river, rpocriverID,'time','river_time');
netcdf.putAtt(nc_river,rpocriverID,'field','river_rpoc, scalar, series')

rponriverID = netcdf.defVar(nc_river,'river_RPON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,rponriverID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
netcdf.putAtt(nc_river,rponriverID,'units','mg N L-1');
netcdf.putAtt(nc_river, rponriverID,'time','river_time');
netcdf.putAtt(nc_river,rponriverID,'field','river_rpon, scalar, series')

rpopriverID = netcdf.defVar(nc_river,'river_RPOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,rpopriverID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
netcdf.putAtt(nc_river,rpopriverID,'units','mg P L-1');
netcdf.putAtt(nc_river, rpopriverID,'time','river_time');
netcdf.putAtt(nc_river,rpopriverID,'field','river_rpop, scalar, series')

salriverID = netcdf.defVar(nc_river,'river_SAL','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,salriverID,'long_name','SALINITY');
netcdf.putAtt(nc_river,salriverID,'units','psu');
netcdf.putAtt(nc_river, salriverID,'time','river_time');
netcdf.putAtt(nc_river,salriverID,'field','river_sal, scalar, series')

sitriverID = netcdf.defVar(nc_river,'river_SIT','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_river,sitriverID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
netcdf.putAtt(nc_river,sitriverID,'units','mg SI L-1');
netcdf.putAtt(nc_river, sitriverID,'time','river_time');
netcdf.putAtt(nc_river,sitriverID,'field','river_sit, scalar, series')

netcdf.close(nc_river)
end