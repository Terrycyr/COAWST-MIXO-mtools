clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check following variables when setting up the script
% 1. Grid name.
% 2. Year
% 3. nontidal_dataset_dir and nontidal_dataset_source
% 4. Output file name
% 5. nfiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%grd_name =  '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
grd_name =  '../Model_grid/ROMS_WFS_Piney.nc';
lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd_name,'Cs_r'));
nfiles = 6;

year = 2021;
%nontidal_dataset_dir = '../Non_tidal_component_preprocessing/HYCOM/Detided_GOM_HYCOM/';
%nontidal_dataset_source = 'GOMu0.04/expt_90.1m000';
nontidal_dataset_dir = '../Non_tidal_component_preprocessing/HYCOM/';
nontidal_dataset_source = 'GLBy0.08/expt_93.0';

load(strcat(nontidal_dataset_dir,'non_tidal_clm_',num2str(year),'.mat'));
origin_date = date_clm + datenum(year,1,1);
ldate = round(length(origin_date)/nfiles);

for i=1:nfiles
    eval(strcat('fn{',num2str(i),'} = ''WFS_',num2str(year),'_clm',num2str(i),'.nc'';'));

    if(i==1)
        if((origin_date(1)-datenum(year,1,1))>0)
            out_date{i} = [datenum(year,1,1);origin_date(1:ldate)];
            date_clm = [0;date_clm];
        else
            out_date{i} = origin_date(1:ldate);
        end
    elseif(i==nfiles)
        if((origin_date(end)-datenum(year+1,1,1))<0)
            out_date{i} = [origin_date(((i-1)*ldate+1):end);datenum(year+1,1,1)];
            date_clm = [date_clm;datenum(year+1,1,1)-datenum(year,1,1)];
        else
            out_date{i} = origin_date(((i-1)*ldate+1):end);
        end
    else
        out_date{i} = origin_date([1:ldate]+(i-1)*ldate);
    end

    for j=1:length(out_date{i})
        delta = abs(origin_date(:)-out_date{i}(j));
        out_i{i}(j) = find(delta==min(delta));
    end
    date_clm_sub{i} = out_date{i}-datenum(year,1,1);

    t_clim(i) = length(date_clm_sub{i});

    create_roms_netcdf_clm_mwUL(fn{i},gn,t_clim(i),nontidal_dataset_source);

    u_clm_sub = u_clm_out(:,:,:,out_i{i});
    v_clm_sub = v_clm_out(:,:,:,out_i{i});
    u2d_clm_sub = u2d_clm_out(:,:,out_i{i});
    v2d_clm_sub = v2d_clm_out(:,:,out_i{i});
    s_clm_sub = s_clm_out(:,:,:,out_i{i});
    temp_clm_sub = temp_clm_out(:,:,:,out_i{i});

    ncwrite(fn{i},'ocean_time',date_clm_sub{i});
    ncwrite(fn{i},'v2d_time',date_clm_sub{i});
    ncwrite(fn{i},'v3d_time',date_clm_sub{i});
    ncwrite(fn{i},'salt_time',date_clm_sub{i});
    ncwrite(fn{i},'temp_time',date_clm_sub{i});
    ncwrite(fn{i},'lon_rho',lon);
    ncwrite(fn{i},'lat_rho',lat);
    ncwrite(fn{i},'u',u_clm_sub);
    ncwrite(fn{i},'v',v_clm_sub);
    ncwrite(fn{i},'ubar',u2d_clm_sub);
    ncwrite(fn{i},'vbar',v2d_clm_sub);
    ncwrite(fn{i},'salt',s_clm_sub);
    ncwrite(fn{i},'temp',temp_clm_sub);
end
