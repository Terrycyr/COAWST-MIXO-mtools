clear all;close all;
%Reading and pre-processing CFSR Variables, and save to mat files.
%Developed by YUren Chen, 2021/03/02, Guangzhou
%Modified for adaption to CFSR and CFSv2, 2021/03/04, Guangzhou
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify Variables
grd = '../Model_grid/ROMS_WFS_new.nc'; %ROMS grid
cfsr_dir = {'../CFSR/1998_2003/wind/','../CFSR/1998_2003/wind/','../CFSR/1998_2003/slp/'}; %CFSR directory for each variables
cfsr_var = {'U_GRD_L103','V_GRD_L103','PRMSL_L101'}; %CFSR variable names
cfsr_deltat = 6; %time interval of the CFSR data, unit: hour
out_var = {'u_out','v_out','slp_out'}; %name of the output variables
origin_date = datenum(2002,1,1,6,0,0):cfsr_deltat/24:datenum(2002,12,31,24,0,0); %origin_date should Cover the Initial time
out_date = datenum(2002,1,1,0,0,0):cfsr_deltat/24:datenum(2002,12,31,24,0,0); %output time range
rot_flag = [1,2,0]; %rotation flag, u==1, v==2, others == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------NO NEED TO CHANGE BELOW--------------------------
year = datevec(out_date(1));
year = year(1);
outname = strcat('atm_forcing_',num2str(year),'.mat');
lon = ncread(grd,'lon_rho');
lat = ncread(grd, 'lat_rho');
angle = ncread(grd,'angle');
year = str2double(datestr(out_date(1), 'yyyy'));

%Read data
read_dat = cell(1,length(out_var));
lon_cfsr = cell(1,length(out_var));
lat_cfsr = cell(1,length(out_var));
for time = origin_date(1):cfsr_deltat/24:origin_date(end) 
    for var_n = 1:length(out_var)
        [cfsr_fname,time_id] = filename_tid(cfsr_dir{var_n},time,cfsr_var{var_n}); 
        tmp = ncread(cfsr_fname,cfsr_var{var_n},[1,1,time_id],[Inf Inf 1]);
        if(time == origin_date(1))
            x = double(ncread(cfsr_fname,'lon'));
            x(x>180) = x(x>180)-360;
            y = double(ncread(cfsr_fname,'lat')); 
            [lat_cfsr{var_n},lon_cfsr{var_n}] = meshgrid(y,x); 
            read_dat{var_n} = tmp;
        else
            tsize_tmp = size(tmp,3);
            tsize_dat = size(read_dat{var_n},3);
            read_dat{var_n}(:,:,[1:tsize_tmp]+tsize_dat) = tmp;
        end
        
    end
end

%Pre-processing for Saving
tinterp_dat = cell(1,length(out_var));
sinterp_dat = cell(1,length(out_var));
for var_n = 1:length(out_var)
    %Temporal interpolation
    tinterp_dat{var_n} = tv3dinterpt(read_dat{var_n}, origin_date, out_date, 'linear');    
    %Spatial interpolation
    input_x = lon_cfsr{var_n};
    input_y = lat_cfsr{var_n};
    out_x = lon;
    out_y = lat;
    sinterp_dat{var_n} = tv3dinterps(tinterp_dat{var_n}, input_x, input_y, out_x, out_y,'linear');   
end

clear read_dat

%Rotation
if(sum(rot_flag) == 3)
    u = sinterp_dat{rot_flag==1};
    v = sinterp_dat{rot_flag==2};
    [sinterp_dat{rot_flag==1}, sinterp_dat{rot_flag==2}] = tv3drot(u, v, angle);
end

%Saving
for var_n = 1:length(out_var)
    eval(strcat(out_var{var_n},'= sinterp_dat{var_n};'));
end

if(~exist(outname,'file'))
    save(outname,"out_date",'-v7.3')
end
for var_n = 1:length(out_var)
    eval(strcat('save(outname,','out_var{var_n},''-append'');'));
end