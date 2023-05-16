clear all; close all;
addpath(path,'../Tidal_component_preprocessing/OTPS2frc/');
addpath(path,'../Tidal_component_preprocessing/OTPS2frc/TMD/');
addpath(path,'../Tidal_component_preprocessing/OTPS2frc/TMD/FUNCTIONS');
addpath(path,'../Tidal_component_preprocessing/OTPS2frc/t_tide_v1.4beta/');

year = 2010;
gfile='ROMS_WFS_10river_grid_v11.nc';
base_date=datenum(year,1,1);
pred_date=(datenum(year+1,1,1)-base_date)/2;
ofile=['WFS_',num2str(year),'_tides_2.nc'];
model_file='../Tidal_component_preprocessing/OTPS2frc/DATA/Model_atlas_v1';
%model_file='../Tidal_component_preprocessing/OTPS2frc/DATA/Model_Mex';
otps2frc_v5(gfile,base_date,pred_date,ofile,model_file,'WFS')