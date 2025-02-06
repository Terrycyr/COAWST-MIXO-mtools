clear all;

%Reading and pre-processing WRF results for ROMS, and save to mat files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify Variables
grd = '/home/ychen/EcoHAB/Model_grid/ROMS_WFS_10river_grid_v10.nc'; %ROMS grid
res_dir = {'../Ian_D2/'};  %directory for WRF results
domain_id = {'02'};
res_var = {'U10','V10','PSFC','SST','SWDNB','SWUPB','LWDNB','RAINC+RAINNC'}; %CFSR variable names
res_deltat = 1; %time interval of the WRF results, unit: hour
out_var = {'u_out','v_out','slp_out','t_out','ssr_d_out','ssr_u_out','slr_d_out','rain_out'}; %name of the output variables
origin_date = datenum(2022,9,27,0,0,0):res_deltat/24:datenum(2022,10,2,0,0,0); %origin_date should Cover the Initial time
out_date = datenum(2022,9,27,0,0,0):1/24:datenum(2022,10,2,24,0,0); %output time range
rot_flag = [1,2,0,0,0,0,0,0]; %rotation flag, u==1, v==2, others == 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------NO NEED TO CHANGE BELOW--------------------------
year = datevec(out_date(1));
year = year(1);
outname = strcat('wrf_forcing_',num2str(year),'.mat');
lon = ncread(grd,'lon_rho');
lat = ncread(grd, 'lat_rho');
angle = ncread(grd,'angle');
year = str2double(datestr(out_date(1), 'yyyy'));

%Read data
read_dat = cell(1,length(out_var));
for time = origin_date(1):res_deltat/24:origin_date(end) 
    for var_n = 1:length(out_var)
        if(length(res_dir)>1)
            fn = [res_dir{var_n},'/wrfout_d',domain_id{var_n},'_',datestr(time,'yyyy-mm-dd_HH:MM:SS')];
        else
            fn = [res_dir{1},'/wrfout_d',domain_id{1},'_',datestr(time,'yyyy-mm-dd_HH:MM:SS')];
        end

        if(~isempty(strfind(res_var{var_n},'+')))
            tag_pos = strfind(res_var{var_n},'+');
            for j = 1:length(tag_pos)
                if(j==1)
                    var_name = res_var{var_n}(1:tag_pos(1)-1);
                    tmp = double(ncread(fn,var_name));
                elseif(j==length(tag_pos))
                    var_name = res_var{var_n}(tag_pos(end)+1:end);
                    tmp = tmp+double(ncread(fn,var_name));
                else
                    var_name = res_var{var_n}(tag_pos(j)+1:tag_pos(j+1)-1);
                    tmp = tmp+double(ncread(fn,var_name));
                end
            end
        else
            tmp = double(ncread(fn,res_var{var_n}));
        end

        if(time == origin_date(1))
            lon_wrf = double(ncread(fn,'XLONG'));
            lon_wrf(lon_wrf>180) = lon_wrf(lon_wrf>180)-360;
            lat_wrf = double(ncread(fn,'XLAT')); 
            read_dat{var_n} = tmp;
        else
            tsize_dat = size(read_dat{var_n},3);
            read_dat{var_n}(:,:,1+tsize_dat) = tmp;
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
    input_x = lon_wrf;
    input_y = lat_wrf;
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