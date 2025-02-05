% script create_roms_rivers
%
% Create a netcdf file that contains river forcing data for ROMS.
% Forcing data consists of:
% 'river_Xposition'  -   'river runoff  XI-positions at RHO-points'
% 'river_Eposition'  -   'river runoff ETA-positions at RHO-points'
% 'river_direction'  -   'river runoff direction'
% 'river_Vshape'     -   'river runoff mass transport vertical profile'
% 'river_transport'  -   'river runoff mass transport'
% 'river_temp'       -   'river runoff potential temperature'
% 'river_salt'       -   'river runoff salinity'
% 'river_mud_'       -   'river runoff suspended sediment concentration'
% 'river_sand_'      -   'river runoff suspended sediment concentration'
%
%
% This m file is set to force rivers for LAKE_SIGNELL for ocean2.2 release.
% Users can adapt this file to their own application.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Begin user input section                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%!!!!!!!!!!!!!!!!!!!!!!!!!Check the position of the rivers!!!!!!!!!!!!!!!!!
%0) Preprocessing of river data.
    addpath(path,'C:\Users\cheny\Desktop\COAWST\Tools\mfiles\mtools');
    
    out_date = 0:1/24:365;

    river_flow = zeros(1,length(out_date));
    river_salt = zeros(1,21,length(out_date));
    river_temp = zeros(1,21,length(out_date));

    river_LDOC = zeros(1,21,length(out_date));
    river_LDON = zeros(1,21,length(out_date));
    river_LDOP = zeros(1,21,length(out_date));
    river_RDOC = zeros(1,21,length(out_date));
    river_RDON = zeros(1,21,length(out_date));
    river_RDOP = zeros(1,21,length(out_date));
    river_NH4T = zeros(1,21,length(out_date));
    river_NO23 = zeros(1,21,length(out_date));
    river_PO4T = zeros(1,21,length(out_date));
    river_LPOC = zeros(1,21,length(out_date));
    river_LPON = zeros(1,21,length(out_date));
    river_LPOP = zeros(1,21,length(out_date));
    river_RPOC = zeros(1,21,length(out_date));
    river_RPON = zeros(1,21,length(out_date));
    river_RPOP = zeros(1,21,length(out_date));
    river_SAL = zeros(1,21,length(out_date));
    river_SIT = zeros(1,21,length(out_date));
    river_DO =  zeros(1,21,length(out_date));

    river_flag = abs(river_flow(end,:))>0;
    river_NH4T(end,:,river_flag) = 0;

%1) Enter name of netcdf forcing file to be created.
%   If it already exists it will be overwritten!!.
    forc_file='bio1dtest_river_bio.nc';

%2) Enter times of river forcings data, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
%     river_time= out_date - out_date(1);      % 3 days -- NEEDS to be in your time code.
     river_time= 0:1/24:365;
%
     num_river_times=length(river_time);     % do not change this.

%3) Enter number of vertical sigma levels in model.
%   This will be same value as entered in mod_param.F
     N=21;

%4) Enter the values of theta_s, theta_b, and Tcline from your *.in file.
%    theta_s = 8;  
%    theta_b = 4;  
%    Tcline =  0; 

%5) Enter value of h, Lm, and Mm.
%   This info can come from a grid file or user supplied here.
%   
%   Are you entering a grid file name (1 = yes, 0 = no)? 
    get_grid = 1;    %<--- put a 1 or 0 here
  
    if (get_grid)
      grid_file = '../Model_grid/ROMS_WFS_Piney.nc';    %<-enter name of grid here
%
% Get some grid info, do not change this.
% 
      netcdf_load(grid_file);
      [LP,MP]=size(h);
      Lm=LP-2;
      Mm=MP-2;
%
    else
      Lm=100;       %<--- else put size of grid here, from mod_param.F
      Mm=20;        %<--- else put size of grid here, from mod_param.F
      LP = Lm+2;    %don't change this.
      MP = Mm+2;    %don't change this.

      % enter depth, same as in ana_grid
      for j=1:MP
        for i=1:LP
          h(j,i)=18-16*(Mm-j)/(Mm-1);
        end
      end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc some grid stuff here - do not change this.
% You go on to step 6.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   L  = Lm+1;
   M  = Mm+1;
   xi_psi  = L;
   xi_rho  = LP;
   xi_u    = L;
   xi_v    = LP;
   eta_psi = M;
   eta_rho = MP;
   eta_u   = MP;
   eta_v   = M;
%
% Don't change this either.  This is from set_scoord.F
% This info is calculated so you know the vertical spacings
% of the grid at startup.
% Go to step 6 now.
%
   hmin=0.1;
   hc=min([hmin,Tcline]);
   sc_r = s_rho;
   sc_w = s_w;
   if(~exist('Cs_r','var')||~exist('Cs_w','var'))
       
       if (theta_s~=0.0)
         cff1=1.0/sinh(theta_s);
         cff2=0.5/tanh(0.5*theta_s);
       end
       sc_w(1)=-1.0;
       Cs_w(1)=-1.0;
       cff=1.0/N;
       for k=1:N
         sc_w(k+1)=cff*(k-N);
         sc_r(k)=cff*((k-N)-0.5);
         if (theta_s~=0)
           Cs_w(k+1)=(1.0-theta_b)*cff1*sinh(theta_s*sc_w(k+1))+   ...
                          theta_b*(cff2*tanh(theta_s*(sc_w(k+1)+0.5))-0.5);
           Cs_r(k)  =(1.0-theta_b)*cff1*sinh(theta_s*sc_r(k))+   ...
                          theta_b*(cff2*tanh(theta_s*(sc_r(k)+0.5))-0.5);
         else
           Cs_w(k+1)=sc_w(k+1);
           Cs_r(k)=sc_r(k);
         end
       end
   end
%
% Assume zeta starts at 0.
%
    for j=1:eta_rho
      for i=1:xi_rho
        zeta(1,j,i) = 0;
      end
    end

%6) Enter number of rivers.
    num_rivers=size(river_flow,1);

%7) Initialize river location and direction parameters.
%   The offline mass sources are located at water grid points. 
%   Use mask_rho to check the position
    lon2 = lon_rho(2:end-1,2:end-1);
    lat2 = lat_rho(2:end-1,2:end-1);
    mask2 = mask_rho(2:end-1,2:end-1);

    river_Xposition=[3]-1;
    river_Eposition=[3]-1;  % num_rivers
    river_direction=[0]; 

    %Check plots
    figure(1);
    hp = pcolor(mask2');  %Check river_direction
    set(hp,'linestyle','none');
    hold on; 
    scatter(river_Xposition,river_Eposition,'g','filled');

%8) Initialize river shape.
    for i=1:num_rivers
      for k=1:N
        %river_Vshape(i,k)=1/N;
        river_Vshape(i,k) = abs(Cs_w(k)-Cs_w(k+1));
      end
      if((sum(river_Vshape(i,:))-1)>0.001)
          'Sum is not 1. Check Vashape!'
          pause;
      end
    end

%9) Initialize river flow.
    river_transport=river_flow; %read in from file
%     river_transport=ones(num_rivers,num_river_times);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END of USER INPUT                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create init file
nc_forc=netcdf.create(forc_file, 'clobber');
if isempty(nc_forc)
	disp([' ## Unable to create ROMS Rivers NetCDF file.'])
	return
end
 
%% Global attributes:
disp(' ## Defining Global Attributes...')
netcdf.putAtt(nc_forc,netcdf.getConstant('NC_GLOBAL'),'type','ROMS Rivers forcing file');
netcdf.putAtt(nc_forc,netcdf.getConstant('NC_GLOBAL'),'history',['Created by ', mfilename ', on ', datestr(now)]);
netcdf.putAtt(nc_forc,netcdf.getConstant('NC_GLOBAL'),'title','ROMS Application')

%% Dimensions:
disp(' ## Defining Dimensions...')

xpsidimID = netcdf.defDim(nc_forc,'xpsi',L);
xrhodimID = netcdf.defDim(nc_forc,'xrho',LP);
xudimID   = netcdf.defDim(nc_forc,'xu',L);
xvdimID   = netcdf.defDim(nc_forc,'xv',LP);

epsidimID = netcdf.defDim(nc_forc,'epsi',M);
erhodimID = netcdf.defDim(nc_forc,'erho',MP);
eudimID   = netcdf.defDim(nc_forc,'eu',MP);
evdimID   = netcdf.defDim(nc_forc,'ev',M);

s_rhodimID = netcdf.defDim(nc_forc,'s_rho',N);
s_wdimID = netcdf.defDim(nc_forc,'s_w',N+1);
riverdimID = netcdf.defDim(nc_forc,'river',num_rivers);
river_timedimID = netcdf.defDim(nc_forc,'river_time',num_river_times);
onedimID = netcdf.defDim(nc_forc,'one',1);

%% Variables and attributes:
disp(' ## Defining Dimensions, Variables, and Attributes...')
 
theta_sID = netcdf.defVar(nc_forc,'theta_s','double',onedimID);
netcdf.putAtt(nc_forc,theta_sID,'long_name','S-coordinate surface control parameter');
netcdf.putAtt(nc_forc,theta_sID,'units','1');

theta_bID = netcdf.defVar(nc_forc,'theta_b','double',onedimID);
netcdf.putAtt(nc_forc,theta_bID,'long_name','S-coordinate bottom control parameter');
netcdf.putAtt(nc_forc,theta_bID,'units','1');

tcline_ID = netcdf.defVar(nc_forc,'Tcline','double',onedimID);
netcdf.putAtt(nc_forc,tcline_ID,'long_name','S-coordinate surface/bottom layer width');
netcdf.putAtt(nc_forc,tcline_ID,'units','meter');

hc_ID = netcdf.defVar(nc_forc,'hc','double',onedimID);
netcdf.putAtt(nc_forc,hc_ID,'long_name','S-coordinate parameter, critical depth');
netcdf.putAtt(nc_forc,hc_ID,'units','meter');

Cs_rID = netcdf.defVar(nc_forc,'Cs_r','double',s_rhodimID);
netcdf.putAtt(nc_forc,Cs_rID,'long_name','S-coordinate stretching curves at RHO-points');
netcdf.putAtt(nc_forc,Cs_rID,'units','1');
netcdf.putAtt(nc_forc,Cs_rID,'valid_min',-1);
netcdf.putAtt(nc_forc,Cs_rID,'valid_max',0);
netcdf.putAtt(nc_forc,Cs_rID,'field','Cs_r, scalar');

Cs_wID = netcdf.defVar(nc_forc,'Cs_w','double',s_wdimID);
netcdf.putAtt(nc_forc,Cs_wID,'long_name','S-coordinate stretching curves at W-points');
netcdf.putAtt(nc_forc,Cs_wID,'units','1');
netcdf.putAtt(nc_forc,Cs_wID,'valid_min',-1);
netcdf.putAtt(nc_forc,Cs_wID,'valid_max',0);
netcdf.putAtt(nc_forc,Cs_wID,'field','Cs_w, scalar');

sc_rID = netcdf.defVar(nc_forc,'sc_r','double',s_rhodimID);
netcdf.putAtt(nc_forc,sc_rID,'long_name','S-coordinate at RHO-points');
netcdf.putAtt(nc_forc,sc_rID,'units','1');
netcdf.putAtt(nc_forc,sc_rID,'valid_min',-1);
netcdf.putAtt(nc_forc,sc_rID,'valid_max',0);
netcdf.putAtt(nc_forc,sc_rID,'field','sc_r, scalar');

sc_wID = netcdf.defVar(nc_forc,'sc_w','double',s_wdimID);
netcdf.putAtt(nc_forc,sc_wID,'long_name','S-coordinate at W-points');
netcdf.putAtt(nc_forc,sc_wID,'units','1');
netcdf.putAtt(nc_forc,sc_wID,'valid_min',-1);
netcdf.putAtt(nc_forc,sc_wID,'valid_max',0);
netcdf.putAtt(nc_forc,sc_wID,'field','sc_w, scalar');

river_ID = netcdf.defVar(nc_forc,'river','double',riverdimID);
netcdf.putAtt(nc_forc,river_ID,'long_name','river_runoff identification number');
netcdf.putAtt(nc_forc,river_ID,'units','nondimensional');
netcdf.putAtt(nc_forc,river_ID,'field','num_rivers, scalar, series');

river_timeID = netcdf.defVar(nc_forc,'river_time','double',river_timedimID);
netcdf.putAtt(nc_forc,river_timeID,'long_name','river time');
netcdf.putAtt(nc_forc,river_timeID,'units','days');
netcdf.putAtt(nc_forc,river_timeID,'field','river_time, scalar, series');

river_XpositionID = netcdf.defVar(nc_forc,'river_Xposition','double',riverdimID);
netcdf.putAtt(nc_forc,river_XpositionID,'long_name','river runoff  XI-positions at RHO-points');
netcdf.putAtt(nc_forc,river_XpositionID,'units','scalar');
netcdf.putAtt(nc_forc,river_XpositionID,'time','river_time');
netcdf.putAtt(nc_forc,river_XpositionID,'field','river runoff XI position, scalar, series');

river_EpositionID = netcdf.defVar(nc_forc,'river_Eposition','double',riverdimID);
netcdf.putAtt(nc_forc,river_EpositionID,'long_name','river runoff  ETA-positions at RHO-points');
netcdf.putAtt(nc_forc,river_EpositionID,'units','scalar');
netcdf.putAtt(nc_forc,river_EpositionID,'time','river_time');
netcdf.putAtt(nc_forc,river_EpositionID,'field','river runoff ETA position, scalar, series');

river_directionID = netcdf.defVar(nc_forc,'river_direction','double',riverdimID);
netcdf.putAtt(nc_forc,river_directionID,'long_name','river runoff direction, XI=0, ETA>0');
netcdf.putAtt(nc_forc,river_directionID,'units','scalar');
netcdf.putAtt(nc_forc,river_directionID,'time','river_time');
netcdf.putAtt(nc_forc,river_directionID,'field','river runoff direction, scalar, series');

river_VshapeID = netcdf.defVar(nc_forc,'river_Vshape','double',[riverdimID s_rhodimID]);
netcdf.putAtt(nc_forc,river_VshapeID,'long_name','river runoff mass transport vertical profile');
netcdf.putAtt(nc_forc,river_VshapeID,'units','scalar');
netcdf.putAtt(nc_forc,river_VshapeID,'time','river_time');
netcdf.putAtt(nc_forc,river_VshapeID,'field','river runoff vertical profile, scalar, series');

river_transportID = netcdf.defVar(nc_forc,'river_transport','double',[riverdimID river_timedimID]);
netcdf.putAtt(nc_forc,river_transportID,'long_name','river runoff mass transport');
netcdf.putAtt(nc_forc,river_transportID,'units','meter^3/s');
netcdf.putAtt(nc_forc,river_transportID,'time','river_time');
netcdf.putAtt(nc_forc,river_transportID,'field','river runoff mass transport, scalar, series');

acriverID = netcdf.defVar(nc_forc,'river_A_C','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,acriverID,'long_name','MX model prey');
netcdf.putAtt(nc_forc,acriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, acriverID,'time','river_time');
netcdf.putAtt(nc_forc,acriverID,'field','river_A_C, scalar, series');

achlcriverID = netcdf.defVar(nc_forc,'river_A_CHLC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,achlcriverID,'long_name','MX model CHLC/C of Prey');
netcdf.putAtt(nc_forc,achlcriverID,'units','gChl/gC');
netcdf.putAtt(nc_forc, achlcriverID,'time','river_time');
netcdf.putAtt(nc_forc,achlcriverID,'field','river_A_CHLC, scalar, series');

ancriverID = netcdf.defVar(nc_forc,'river_A_NC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,ancriverID,'long_name','MX model N:C of Prey');
netcdf.putAtt(nc_forc,ancriverID,'units','gN/gC');
netcdf.putAtt(nc_forc, ancriverID,'time','river_time');
netcdf.putAtt(nc_forc,ancriverID,'field','river_A_NC, scalar, series');

apcriverID = netcdf.defVar(nc_forc,'river_A_PC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,apcriverID,'long_name','MX model P:C of Prey');
netcdf.putAtt(nc_forc,apcriverID,'units','gP/gC');
netcdf.putAtt(nc_forc, apcriverID,'time','river_time');
netcdf.putAtt(nc_forc,apcriverID,'field','river_A_PC, scalar, series');

bsiriverID = netcdf.defVar(nc_forc,'river_BSI','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,bsiriverID,'long_name','BIOGENIC SILICA');
netcdf.putAtt(nc_forc,bsiriverID,'units','mg SI L-1');
netcdf.putAtt(nc_forc, bsiriverID,'time','river_time');
netcdf.putAtt(nc_forc,bsiriverID,'field','river_BSI, scalar, series');

doriverID = netcdf.defVar(nc_forc,'river_DO','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,doriverID,'long_name','DISSOLVED OXYGEN');
netcdf.putAtt(nc_forc,doriverID,'units','mg O2 L-1');
netcdf.putAtt(nc_forc, doriverID,'time','river_time');
netcdf.putAtt(nc_forc,doriverID,'field','river_DO, scalar, series');

exdocriverID = netcdf.defVar(nc_forc,'river_EXDOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,exdocriverID,'long_name','ALGAL EXUDATE - DISSOLVED ORGANIC CARBON');
netcdf.putAtt(nc_forc,exdocriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, exdocriverID,'time','river_time');
netcdf.putAtt(nc_forc,exdocriverID,'field','river_exdoc, scalar, series');

ldocriverID = netcdf.defVar(nc_forc,'river_LDOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,ldocriverID,'long_name','LABILE DISSOLVED ORGANIC CARBON (LDOC)');
netcdf.putAtt(nc_forc,ldocriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, ldocriverID,'time','river_time');
netcdf.putAtt(nc_forc,ldocriverID,'field','river_ldoc, scalar, series');

ldonriverID = netcdf.defVar(nc_forc,'river_LDON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,ldonriverID,'long_name','LABILE DISSOLVED ORGANIC NITROGEN (LDON)');
netcdf.putAtt(nc_forc,ldonriverID,'units','mg N L-1');
netcdf.putAtt(nc_forc, ldonriverID,'time','river_time');
netcdf.putAtt(nc_forc,ldonriverID,'field','river_ldon, scalar, series');

ldopriverID = netcdf.defVar(nc_forc,'river_LDOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,ldopriverID,'long_name','LABILE DISSOLVED ORGANIC PHOSPHORUS (LDOP)');
netcdf.putAtt(nc_forc,ldopriverID,'units','mg P L-1');
netcdf.putAtt(nc_forc, ldopriverID,'time','river_time');
netcdf.putAtt(nc_forc,ldopriverID,'field','river_ldop, scalar, series');

lpocriverID = netcdf.defVar(nc_forc,'river_LPOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,lpocriverID,'long_name','LABILE PARTICULATE ORGANIC CARBON (LPOC)');
netcdf.putAtt(nc_forc,lpocriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, lpocriverID,'time','river_time');
netcdf.putAtt(nc_forc,lpocriverID,'field','river_lpoc, scalar, series');

lponriverID = netcdf.defVar(nc_forc,'river_LPON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,lponriverID,'long_name','LABILE PARTICULATE ORGANIC NITROGEN (LPON)');
netcdf.putAtt(nc_forc,lponriverID,'units','mg N L-1');
netcdf.putAtt(nc_forc, lponriverID,'time','river_time');
netcdf.putAtt(nc_forc,lponriverID,'field','river_lpon, scalar, series');

lpopriverID = netcdf.defVar(nc_forc,'river_LPOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,lpopriverID,'long_name','LABILE PARTICULATE ORGANIC PHOSPHORUS (LPOP)');
netcdf.putAtt(nc_forc,lpopriverID,'units','mg P L-1');
netcdf.putAtt(nc_forc, lpopriverID,'time','river_time');
netcdf.putAtt(nc_forc,lpopriverID,'field','river_lpop, scalar, series');

mavguriverID = netcdf.defVar(nc_forc,'river_M_AVGU','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mavguriverID,'long_name','MX model mixotroph average growth rate');
netcdf.putAtt(nc_forc,mavguriverID,'units','C/C/d');
netcdf.putAtt(nc_forc, mavguriverID,'time','river_time');
netcdf.putAtt(nc_forc,mavguriverID,'field','river_mavgu, scalar, series');

mcriverID = netcdf.defVar(nc_forc,'river_M_C','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mcriverID,'long_name','MX model core mixotroph body biomass carbon');
netcdf.putAtt(nc_forc,mcriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, mcriverID,'time','river_time');
netcdf.putAtt(nc_forc,mcriverID,'field','river_mc, scalar, series');

mchlcriverID = netcdf.defVar(nc_forc,'river_M_CHLC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mchlcriverID,'long_name','MX model mixotroph core Chla:C');
netcdf.putAtt(nc_forc,mchlcriverID,'units','gChl/gC');
netcdf.putAtt(nc_forc, mchlcriverID,'time','river_time');
netcdf.putAtt(nc_forc,mchlcriverID,'field','river_mchlc, scalar, series');

mfcriverID = netcdf.defVar(nc_forc,'river_M_FC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mfcriverID,'long_name','MX model mixotroph gut relative to body biomass');
netcdf.putAtt(nc_forc,mfcriverID,'units','gC/gC');
netcdf.putAtt(nc_forc, mfcriverID,'time','river_time');
netcdf.putAtt(nc_forc,mfcriverID,'field','river_mfc, scalar, series');

mfchlcriverID = netcdf.defVar(nc_forc,'river_M_FCHLC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mfchlcriverID,'long_name','MX model (ChlA in the gut)/(the core mixotroph biomass)');
netcdf.putAtt(nc_forc,mfchlcriverID,'units','gChl in the gut /gC');
netcdf.putAtt(nc_forc, mfchlcriverID,'time','river_time');
netcdf.putAtt(nc_forc,mfchlcriverID,'field','river_mfchlc, scalar, series');

mncriverID = netcdf.defVar(nc_forc,'river_M_NC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mncriverID,'long_name','MX model mixotroph N:C');
netcdf.putAtt(nc_forc,mncriverID,'units','gN/gC');
netcdf.putAtt(nc_forc, mncriverID,'time','river_time');
netcdf.putAtt(nc_forc,mncriverID,'field','river_mnc, scalar, series');

mpcriverID = netcdf.defVar(nc_forc,'river_M_PC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mpcriverID,'long_name','MX model mixotroph P:C');
netcdf.putAtt(nc_forc,mpcriverID,'units','gP/gC');
netcdf.putAtt(nc_forc, mpcriverID,'time','river_time');
netcdf.putAtt(nc_forc,mpcriverID,'field','river_mpc, scalar, series');

mforc_filecriverID = netcdf.defVar(nc_forc,'river_M_FNC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mforc_filecriverID,'long_name','MX model N:C in mixotroph gut (compared to M_C)');
netcdf.putAtt(nc_forc,mforc_filecriverID,'units','gN/gC');
netcdf.putAtt(nc_forc, mforc_filecriverID,'time','river_time');
netcdf.putAtt(nc_forc,mforc_filecriverID,'field','river_mforc_filec, scalar, series');

mfpcriverID = netcdf.defVar(nc_forc,'river_M_FPC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,mfpcriverID,'long_name','MX model P:C in mixotroph gut (compared to M_C)');
netcdf.putAtt(nc_forc,mfpcriverID,'units','gP/gC');
netcdf.putAtt(nc_forc, mfpcriverID,'time','river_time');
netcdf.putAtt(nc_forc,mfpcriverID,'field','river_mfpc, scalar, series');

nh4triverID = netcdf.defVar(nc_forc,'river_NH4T','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,nh4triverID,'long_name','TOTAL AMMONIA  = Algal N plus dissolved NH4');
netcdf.putAtt(nc_forc,nh4triverID,'units','mg N L-1');
netcdf.putAtt(nc_forc, nh4triverID,'time','river_time');
netcdf.putAtt(nc_forc,nh4triverID,'field','river_nh4t, scalar, series');

no23riverID = netcdf.defVar(nc_forc,'river_NO23','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,no23riverID,'long_name','NITRITE + NITRATE');
netcdf.putAtt(nc_forc,no23riverID,'units','mg N L-1');
netcdf.putAtt(nc_forc, no23riverID,'time','river_time');
netcdf.putAtt(nc_forc,no23riverID,'field','river_no23, scalar, series')

o2eqriverID = netcdf.defVar(nc_forc,'river_O2EQ','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,o2eqriverID,'long_name','AQUEOUS SOD');
netcdf.putAtt(nc_forc,o2eqriverID,'units','mg O2 L-1');
netcdf.putAtt(nc_forc, o2eqriverID,'time','river_time');
netcdf.putAtt(nc_forc,o2eqriverID,'field','river_o2eq, scalar, series')

phyt1riverID = netcdf.defVar(nc_forc,'river_PHYT1','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,phyt1riverID,'long_name','WINTER DIATOMS (PHYT1)');
netcdf.putAtt(nc_forc,phyt1riverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, phyt1riverID,'time','river_time');
netcdf.putAtt(nc_forc,phyt1riverID,'field','river_phyt1, scalar, series')

phyt2riverID = netcdf.defVar(nc_forc,'river_PHYT2','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,phyt2riverID,'long_name','SUMMER ASSEMBLAGE (PHYT2)');
netcdf.putAtt(nc_forc,phyt2riverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, phyt2riverID,'time','river_time');
netcdf.putAtt(nc_forc,phyt2riverID,'field','river_phyt2, scalar, series')

phyt3riverID = netcdf.defVar(nc_forc,'river_PHYT3','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,phyt3riverID,'long_name','FALL ASSEMBLAGE (PHYT3)');
netcdf.putAtt(nc_forc,phyt3riverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, phyt3riverID,'time','river_time');
netcdf.putAtt(nc_forc,phyt3riverID,'field','river_phyt3, scalar, series')

po4triverID = netcdf.defVar(nc_forc,'river_PO4T','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,po4triverID,'long_name','TOTAL DISSOLVED INORGANIC PHOSPHORUS (PO4T) = Algal P plus dissolved PO4');
netcdf.putAtt(nc_forc,po4triverID,'units','mg P L-1');
netcdf.putAtt(nc_forc, po4triverID,'time','river_time');
netcdf.putAtt(nc_forc,po4triverID,'field','river_po4t, scalar, series')

rdocriverID = netcdf.defVar(nc_forc,'river_RDOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,rdocriverID,'long_name','REFRACTORY DISSOLVED ORGANIC CARBON (RDOC)');
netcdf.putAtt(nc_forc,rdocriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, rdocriverID,'time','river_time');
netcdf.putAtt(nc_forc,rdocriverID,'field','river_rdoc, scalar, series')

rdonriverID = netcdf.defVar(nc_forc,'river_RDON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,rdonriverID,'long_name','REFRACTORY DISSOLVED ORGANIC NITROGEN (RDON)');
netcdf.putAtt(nc_forc,rdonriverID,'units','mg N L-1');
netcdf.putAtt(nc_forc, rdonriverID,'time','river_time');
netcdf.putAtt(nc_forc,rdonriverID,'field','river_rdon, scalar, series')

rdopriverID = netcdf.defVar(nc_forc,'river_RDOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,rdopriverID,'long_name','REFRACTORY DISSOLVED ORGANIC PHOSPHORUS (RDOP)');
netcdf.putAtt(nc_forc,rdopriverID,'units','mg P L-1');
netcdf.putAtt(nc_forc, rdopriverID,'time','river_time');
netcdf.putAtt(nc_forc,rdopriverID,'field','river_rdop, scalar, series')

redocriverID = netcdf.defVar(nc_forc,'river_REDOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,redocriverID,'long_name','REACTIVE DISSOLVED ORGANIC CARBON - CSO/WWTP');
netcdf.putAtt(nc_forc,redocriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, redocriverID,'time','river_time');
netcdf.putAtt(nc_forc,redocriverID,'field','river_redoc, scalar, series')

repocriverID = netcdf.defVar(nc_forc,'river_REPOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,repocriverID,'long_name','REACTIVE PARTICULATE ORGANIC CARBON - CSO/WWTP');
netcdf.putAtt(nc_forc,repocriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, repocriverID,'time','river_time');
netcdf.putAtt(nc_forc,repocriverID,'field','river_repoc, scalar, series')

rpocriverID = netcdf.defVar(nc_forc,'river_RPOC','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,rpocriverID,'long_name','REFRACTORY PARTICULATE ORGANIC CARBON (RPOC)');
netcdf.putAtt(nc_forc,rpocriverID,'units','mg C L-1');
netcdf.putAtt(nc_forc, rpocriverID,'time','river_time');
netcdf.putAtt(nc_forc,rpocriverID,'field','river_rpoc, scalar, series')

rponriverID = netcdf.defVar(nc_forc,'river_RPON','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,rponriverID,'long_name','REFRACTORY PARTICULATE ORGANIC NITROGEN (RPON)');
netcdf.putAtt(nc_forc,rponriverID,'units','mg N L-1');
netcdf.putAtt(nc_forc, rponriverID,'time','river_time');
netcdf.putAtt(nc_forc,rponriverID,'field','river_rpon, scalar, series')

rpopriverID = netcdf.defVar(nc_forc,'river_RPOP','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,rpopriverID,'long_name','REFRACTORY PARTICULATE ORGANIC PHOSPHORUS (RPOP)');
netcdf.putAtt(nc_forc,rpopriverID,'units','mg P L-1');
netcdf.putAtt(nc_forc, rpopriverID,'time','river_time');
netcdf.putAtt(nc_forc,rpopriverID,'field','river_rpop, scalar, series')

salriverID = netcdf.defVar(nc_forc,'river_SAL','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,salriverID,'long_name','SALINITY');
netcdf.putAtt(nc_forc,salriverID,'units','psu');
netcdf.putAtt(nc_forc, salriverID,'time','river_time');
netcdf.putAtt(nc_forc,salriverID,'field','river_sal, scalar, series')

sitriverID = netcdf.defVar(nc_forc,'river_SIT','float',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,sitriverID,'long_name','TOTAL INORGANIC SILICA = Algal SI plus dissolved SI');
netcdf.putAtt(nc_forc,sitriverID,'units','mg SI L-1');
netcdf.putAtt(nc_forc, sitriverID,'time','river_time');
netcdf.putAtt(nc_forc,sitriverID,'field','river_sit, scalar, series')
 
netcdf.close(nc_forc)

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')
const_initialize(forc_file,0.);

ncwrite(forc_file,'theta_s',theta_s);
ncwrite(forc_file,'theta_b',theta_b);
ncwrite(forc_file,'Tcline',Tcline);
ncwrite(forc_file,'Cs_r',Cs_r);
ncwrite(forc_file,'Cs_w',Cs_w);
ncwrite(forc_file,'sc_w',sc_w);
ncwrite(forc_file,'sc_r',sc_r);
ncwrite(forc_file,'hc',hc);
ncwrite(forc_file,'river',[1:num_rivers]);

ncwrite(forc_file,'river_time',river_time);
ncwrite(forc_file,'river_Xposition',river_Xposition);
ncwrite(forc_file,'river_Eposition',river_Eposition);
ncwrite(forc_file,'river_direction',river_direction);
ncwrite(forc_file,'river_Vshape',river_Vshape);
ncwrite(forc_file,'river_transport',river_transport);

ncwrite(forc_file,'river_LDOC',river_LDOC);
ncwrite(forc_file,'river_LDON',river_LDON);
ncwrite(forc_file,'river_LDOP',river_LDOP);
ncwrite(forc_file,'river_RDOC',river_RDOC);
ncwrite(forc_file,'river_RDON',river_RDON);
ncwrite(forc_file,'river_RDOP',river_RDOP);
ncwrite(forc_file,'river_NH4T',river_NH4T);
ncwrite(forc_file,'river_NO23',river_NO23);
ncwrite(forc_file,'river_PO4T',river_PO4T);
ncwrite(forc_file,'river_LPOC',river_LPOC);
ncwrite(forc_file,'river_LPON',river_LPON);
ncwrite(forc_file,'river_LPOP',river_LPOP);
ncwrite(forc_file,'river_RPOC',river_RPOC);
ncwrite(forc_file,'river_RPON',river_RPON);
ncwrite(forc_file,'river_RPOP',river_RPOP);
ncwrite(forc_file,'river_SAL',river_SAL);
ncwrite(forc_file,'river_SIT',river_SIT);
ncwrite(forc_file,'river_DO',river_DO);

%close file
disp(['created ', forc_file])




