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
% In your project.h you should have 
% #define TS_PSOURCE
% #define UV_PSOURCE
% #undef  ANA_PSOURCE
%
% This m file is set to force rivers for LAKE_SIGNELL for ocean2.2 release.
% Users can adapt this file to their own application.
%
% jcw 5-25-2005
% jcw 21March2014  to use native matlab netcdf.
% cyr 3-19-2002  
% cyr 6-18-2002  include plot section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Begin user input section                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%0) Preprocessing of river data.
    addpath(path,'../River_preprocessing/USGS_NWIS/');
    dye_flag = 1;
    out_date = datenum(2002,1,1):datenum(2002,12,31,24,0,0);
    year = datevec(out_date(1));
    year = year(1);
    preprocess_river_nwis;

    if(dye_flag==1)
    load(strcat('../Water_Atlas/WA_river_bnd_',num2str(year),'.mat'));
    [a,b,c] = size(river_NO23);
    for i=1:b
        river_tracer01(:,b-i+1,:) = river_NO23(:,i,:);
    end
%    river_tracer01(1:3,:,:) = 0;
    end

%1) Enter name of netcdf forcing file to be created.
%   If it already exists it will be overwritten!!.
    forc_file='WFS_2002_river.nc';

%2) Enter times of river forcings data, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
     river_time= out_date - out_date(1);      % 3 days -- NEEDS to be in your time code.
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
      grid_file = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';    %<-enter name of grid herec
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
   num_rivers=11;

%7) Initialize river location and direction parameters.
%   Currently, the direction can be along XI-direction
%   (river-direction = 0) or along ETA-direction (river_direction > 0).  The mass sources are
%   located at U- or V-points so the grid locations should range from
%   1 =< river_Xposition =< L  and  1 =< river_Eposition =< M
    lon_u2 = lon_u(:,2:end-1);
    lat_u2 = lat_u(:,2:end-1);
    lon_v2 = lon_v(2:end-1,:);
    lat_v2 = lat_v(2:end-1,:);
    mask2 = mask_rho(2:end-1,2:end-1);
    mask_u2 = mask_u(:,2:end-1);
    mask_v2 = mask_v(2:end-1,:);
   
    %river_Xposition=[51;62;75;132;150;161;181;195;208;252]';
    %river_Eposition=[173;191;190;169;179;160;164;140;178;173]'+20;  % num_rivers
    %river_direction=[1;1;1;1;1;0;1;0;1;1]'; 
    
    river_Xposition=[131 149 176 225 246 260 282 303 314 365 216]-1;
    river_Eposition=[170 177 171 143 156 132 136 111 157 155 120]-1;  % num_rivers
    river_direction=[1;1;1;1;1;0;0;0;1;1;0]'; 
    %Check plots
    figure(1);
    hp = pcolor(mask2');  %Check river_direction
    set(hp,'linestyle','none');
    figure(2);
    hp = pcolor(mask_u2'); % left-upper corner ij
    set(hp,'linestyle','none');
    hold on; 
    scatter(river_Xposition(river_direction==0),river_Eposition(river_direction==0),'g','filled');
    figure(3)
    hp = pcolor(mask_v2'); % left-bottom corner ij
    set(hp,'linestyle','none');
    hold on;
    scatter(river_Xposition(river_direction==1),river_Eposition(river_direction==1),'g','filled');
    
    
    
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
     river_flow([4:7 10:11],:) = river_flow([4:7 10:11],:)*1.2;
     river_transport=river_flow; %read in from file
%      river_transport=ones(num_rivers,num_river_times);

%10) Time series of river temp and salt.
    for time=1:num_river_times
      for k=1:N
        for i=1:num_rivers
            if(~exist('river_temp','var'))
              river_temp(i,k,time)=28.5;
            elseif(size(river_temp,2)~=N)
              t_tmp(i,k,time) = river_temp(i,time);
            end
            if(~exist('river_salt','var'))
              river_salt(i,k,time)=0;
            elseif(size(river_salt,2)~=N)
              s_tmp(i,k,time) = river_salt(i,time);
            end
        end
      end
    end
    if(exist('t_tmp','var'))
        river_temp = t_tmp;
    end
    if(exist('s_tmp','var'))
        river_salt = s_tmp*0.0;
    end
        

%11) Enter number of mud sediments (NCS) and number of sand sediments (NNS).
%   These values should be the same as in mod_param.F
    NCS = 0;   %number of cohesive sed classes
    NNS = 0;   %number of non-cohesive sed classes
%
% calc sed parameters. Do not alter.
%
   NAT= 2;  %assume temp + salt are active
   NST = NCS + NNS;     % total number of sed tracers.
   NT = NAT+NST;        % total number of tracers.

%12) Sediment class properties (in order, mud first then sand).
%  These values should coincide with your sediment.in file.
  mud_Srho=ones(1,NCS)*2650;                      %kg m-3, NCS values
  mud_Sd50=[0.50 0.125 0.03125 0.03125 0.03125]/1000;     %m,      NCS values
  mud_Wsed=[40.0 5.0 0.62 0.62 0.62]/1000;             %m s-1,  NCS values
  mud_tau_ce=[0.50 0.10 0.05 0.05 0.05];           %N m-2,  NCS values
  mud_Erate=[3 3 3 30 30]*1e-4;             %kg m-2 s-1, NCS values
  sand_Srho=ones(1,NNS)*2650;       %kg m-3, NNS values
  sand_Sd50=[1.0]/1000;             %m,      NNS values
  sand_Wsed=[1.0]/1000;             %m s-1,  NNS values
  sand_tau_ce=[0.07];               %N m-2,  NNS values
  sand_Erate=[1]*1e-5;              %kg m-2 s-1, NNS values
%
% make some combined arrays.  Do not alter.
%
  Srho=  [mud_Srho,sand_Srho];
  Sd50=  [mud_Sd50,sand_Sd50];
  Wsed=  [mud_Wsed,sand_Wsed];
  tau_ce=[mud_tau_ce,sand_tau_ce];
  Erate= [mud_Erate,sand_Erate];


%13) Time series of river mud and sand.
%
% mud.
  display('Initializing river sediments.')
%
  for idmud=1:NCS
    count=['0',num2str(idmud)];
    count=count(end-1:end);
    for time=1:num_river_times
      for k=1:N
        for i=1:num_rivers
          eval(['river_mud_',count,'(i,k,time) = 0;'])               %mud conc in water column
        end
      end
    end
  end
%
% sand.
%
  for isand=1:NNS
    count=['0',num2str(isand)];
    count=count(end-1:end);
    for time=1:num_river_times
      for k=1:N
        for i=1:num_rivers
          eval(['river_sand_',count,'(i,k,time) = 0;'])               %sand conc in water column
        end
      end
    end
  end

%14) plotting

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
tracerdimID = netcdf.defDim(nc_forc,'tracet',NT);
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

river_tempID = netcdf.defVar(nc_forc,'river_temp','double',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,river_tempID,'long_name','river runoff potential temperature');
netcdf.putAtt(nc_forc,river_tempID,'units','Celsius');
netcdf.putAtt(nc_forc,river_tempID,'time','river_time');
netcdf.putAtt(nc_forc,river_tempID,'field','river temperature, scalar, series');

river_saltID = netcdf.defVar(nc_forc,'river_salt','double',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,river_saltID,'long_name','river runoff salinity');
netcdf.putAtt(nc_forc,river_saltID,'units','PSU');
netcdf.putAtt(nc_forc,river_saltID,'time','river_time');
netcdf.putAtt(nc_forc,river_saltID,'field','river salinity, scalar, series');

if(dye_flag==1)
river_dye01ID = netcdf.defVar(nc_forc,'river_dye_01','double',[riverdimID s_rhodimID river_timedimID]);
netcdf.putAtt(nc_forc,river_dye01ID,'long_name','river runoff tracer 01');
netcdf.putAtt(nc_forc,river_dye01ID,'units','mg/L');
netcdf.putAtt(nc_forc,river_dye01ID,'time','river_time');
netcdf.putAtt(nc_forc,river_dye01ID,'field','river tracer 01, scalar, series');
end

for mm=1:NCS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['river_mud_',count,'ID = netcdf.defVar(nc_forc,''river_mud_',count,''',''double'',[riverdimID s_rhodimID river_timedimID]);'])
  eval(['netcdf.putAtt(nc_forc,river_mud_',count,'ID,''long_name'',''river runoff suspended sediment concentration, size class ',count,''');'])
  eval(['netcdf.putAtt(nc_forc,river_mud_',count,'ID,''units'',''kilogram meter-3'');'])
  eval(['netcdf.putAtt(nc_forc,river_mud_',count,'ID,''time'',''river_time'');'])
  eval(['netcdf.putAtt(nc_forc,river_mud_',count,'ID,''field'',''river runoff mud_',count,', scalar, series'');'])    
end
for mm=1:NNS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['river_sand_',count,'ID = netcdf.defVar(nc_forc,''river_sand_',count,''',''double'',[riverdimID s_rhodimID river_timedimID]);'])
  eval(['netcdf.putAtt(nc_forc,river_sand_',count,'ID,''long_name'',''river runoff suspended sediment concentration, size class ',count,''');'])
  eval(['netcdf.putAtt(nc_forc,river_sand_',count,'ID,''units'',''kilogram meter-3'');'])
  eval(['netcdf.putAtt(nc_forc,river_sand_',count,'ID,''time'',''river_time'');'])
  eval(['netcdf.putAtt(nc_forc,river_sand_',count,'ID,''field'',''river runoff sand_',count,', scalar, series'');'])    
end
netcdf.close(nc_forc)

%now write the data from the arrays to the netcdf file
disp(' ## Filling Variables in netcdf file with data...')

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
ncwrite(forc_file,'river_temp',river_temp);
ncwrite(forc_file,'river_salt',river_salt);

if(dye_flag==1)
ncwrite(forc_file,'river_dye_01',river_tracer01);
end

for mm=1:NCS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['ncwrite(forc_file,''river_mud_',count,''',river_mud_',count,');'])    %mud conc in water column
end
for mm=1:NNS
  count=['00',num2str(mm)];
  count=count(end-1:end);
  eval(['ncwrite(forc_file,''river_sand_',count,''',river_sand_',count,');'])  %sand conc in water column
end

%close file
disp(['created ', forc_file])




