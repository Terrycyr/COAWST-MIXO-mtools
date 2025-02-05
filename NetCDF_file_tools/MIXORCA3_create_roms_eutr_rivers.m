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
%
% jcw 5-25-2005
% jcw 21March2014  to use native matlab netcdf.
% cyr 3-19-2021  
% cyr 6-18-2021  include plot section.
% cyr 3-01-2022  MIXORCA version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Begin user input section                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%0) Preprocessing of river data.
    addpath(path,'C:\Users\cheny\Desktop\EcoHAB\NC_file_generation');
    addpath(path,'../River_preprocessing/USGS_NWIS/');
    out_date = datenum(2005,1,1):datenum(2006,12,31,24,0,0);
    time_ref = datenum(2005,6,1,0,0,0);
    year = datevec(out_date(1));
    year = year(1);
    preprocess_river_nwis;
    wa1 = load(strcat('../Water_Atlas/WA_river_bnd_',num2str(year),'.mat'));
    wa2 = load(strcat('../Water_Atlas/WA_river_bnd_',num2str(year+1),'.mat'));
    
    river_DO = wa1.river_DO;
    river_SAL = wa1.river_SAL;
    river_LDOC = wa1.river_LDOC;
    river_LDON = wa1.river_LDON;
    river_LDOP = wa1.river_LDOP;
    river_LPOC = wa1.river_LPOC;
    river_LPON = wa1.river_LPON;
    river_LPOP = wa1.river_LPOP;
    river_RDOC = wa1.river_RDOC;
    river_RDON = wa1.river_RDON;
    river_RDOP = wa1.river_RDOP;
    river_RPOC = wa1.river_RPOC;
    river_RPON = wa1.river_RPON;
    river_RPOP = wa1.river_RPOP;
    river_NH4T = wa1.river_NH4T;
    river_NO23 = wa1.river_NO23;
    river_PO4T = wa1.river_PO4T;
    river_SIT = wa1.river_SIT;

    pos1 = size(river_DO,3);
    pos2 = pos1+size(wa2.river_DO,3)-1;

    river_DO(:,:,pos1:pos2) = wa2.river_DO;
    river_SAL(:,:,pos1:pos2) = wa2.river_SAL;
    river_LDOC(:,:,pos1:pos2) = wa2.river_LDOC;
    river_LDON(:,:,pos1:pos2) = wa2.river_LDON;
    river_LDOP(:,:,pos1:pos2) = wa2.river_LDOP;
    river_LPOC(:,:,pos1:pos2) = wa2.river_LPOC;
    river_LPON(:,:,pos1:pos2) = wa2.river_LPON;
    river_LPOP(:,:,pos1:pos2) = wa2.river_LPOP;
    river_RDOC(:,:,pos1:pos2) = wa2.river_RDOC;
    river_RDON(:,:,pos1:pos2) = wa2.river_RDON;
    river_RDOP(:,:,pos1:pos2) = wa2.river_RDOP;
    river_RPOC(:,:,pos1:pos2) = wa2.river_RPOC;
    river_RPON(:,:,pos1:pos2) = wa2.river_RPON;
    river_RPOP(:,:,pos1:pos2) = wa2.river_RPOP;
    river_NH4T(:,:,pos1:pos2) = wa2.river_NH4T;
    river_NO23(:,:,pos1:pos2) = wa2.river_NO23;
    river_PO4T(:,:,pos1:pos2) = wa2.river_PO4T;
    river_SIT(:,:,pos1:pos2) = wa2.river_SIT;

    Nutri_adjust = 1.0;

%1) Enter name of netcdf forcing file to be created.
%   If it already exists it will be overwritten!!.
    forc_file='WFS_2005_2006_river_bio_mixo.nc';

%2) Enter times of river forcings data, in seconds.
%   This time needs to be consistent with model time (ie dstart and time_ref).
%   See *.in files for more detail. 
%     river_time= out_date - out_date(1);      % 3 days -- NEEDS to be in your time code.
     river_time= out_date - time_ref;
%
     num_river_times=length(river_time);     % do not change this.

%3) Enter number of vertical sigma levels in model.
%   This will be same value as entered in mod_param.F
     grd_name = '../Model_grid/ROMS_WFS_new.nc';    %<-enter name of grid here
    
     lon_rho = ncread(grd_name,'lon_rho');
     lat_rho = ncread(grd_name,'lat_rho');
     mask_rho = ncread(grd_name,'mask_rho');
     Cs_r = ncread(grd_name,'Cs_r');
     Cs_w = ncread(grd_name,'Cs_w');
     sc_r = ncread(grd_name,'s_rho');
     sc_w = ncread(grd_name,'s_w');
     hc = ncread(grd_name,'hc');
     theta_s = ncread(grd_name,'theta_s');
     theta_b = ncread(grd_name,'theta_b');
     Tcline = ncread(grd_name,'Tcline');
     gn = struct('lon_rho',lon_rho);
     gn.N =length(Cs_r);
     N = gn.N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calc some grid stuff here - do not change this.
% You go on to step 6.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4) Enter number of rivers.
     num_rivers=size(river_flow,1);


%5) Initialize river location and direction parameters.
%   The offline mass sources are located at water grid points. 
%   Use mask_rho to check the position

    lon2 = lon_rho(2:end-1,2:end-1);
    lat2 = lat_rho(2:end-1,2:end-1);
    mask2 = mask_rho(2:end-1,2:end-1);

    river_Xposition=[131 149 176 225 246 259 281 303 314 365 216]-1;
    river_Eposition=[169 176 170 142 155 132 136 111 156 154 120]-1;  % num_rivers
    river_direction=[1;1;1;1;1;0;0;0;1;1;0]'; 

    %Check plots
    figure(1);
    hp = pcolor(mask2');  %Check river_direction
    set(hp,'linestyle','none');
    hold on; 
    scatter(river_Xposition,river_Eposition,'g','filled');

%6) Initialize river shape.
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

%7) Initialize river flow.
    river_transport=river_flow; %read in from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  END of USER INPUT                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create netcdf file
create_roms_netcdf_river(forc_file,gn,num_rivers,num_river_times)

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

ncwrite(forc_file,'river_LDOC',river_LDOC*Nutri_adjust);
ncwrite(forc_file,'river_LDON',river_LDON*Nutri_adjust);
ncwrite(forc_file,'river_LDOP',river_LDOP*Nutri_adjust);
ncwrite(forc_file,'river_RDOC',river_RDOC*Nutri_adjust);
ncwrite(forc_file,'river_RDON',river_RDON*Nutri_adjust);
ncwrite(forc_file,'river_RDOP',river_RDOP*Nutri_adjust);
ncwrite(forc_file,'river_NH4T',river_NH4T*Nutri_adjust);
ncwrite(forc_file,'river_NO23',river_NO23*Nutri_adjust);
ncwrite(forc_file,'river_PO4T',river_PO4T*Nutri_adjust);
ncwrite(forc_file,'river_LPOC',river_LPOC*Nutri_adjust);
ncwrite(forc_file,'river_LPON',river_LPON*Nutri_adjust);
ncwrite(forc_file,'river_LPOP',river_LPOP*Nutri_adjust);
ncwrite(forc_file,'river_RPOC',river_RPOC*Nutri_adjust);
ncwrite(forc_file,'river_RPON',river_RPON*Nutri_adjust);
ncwrite(forc_file,'river_RPOP',river_RPOP*Nutri_adjust);
ncwrite(forc_file,'river_SAL',river_SAL);
ncwrite(forc_file,'river_SIT',river_SIT*Nutri_adjust);
ncwrite(forc_file,'river_DO',river_DO);

%close file
disp(['created ', forc_file])




