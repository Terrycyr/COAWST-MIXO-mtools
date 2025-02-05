clear all;

%Merge WRF results and CFSR data, and save to mat files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify Variables
grd = '../Model_grid/ROMS_WFS_10river_grid_v10.nc'; %ROMS grid
cfsr_dir = 'C:\Users\cheny\Desktop\EcoHAB\Atmosphere_forcing_preprocessing'; %CFSR directory for each variables
merge_var = {'u_out','v_out','slp_out','t_out','ssr_d_out','ssr_u_out','slr_d_out','hum_out','rain_out'}; %name of the output variables
merge_date = datenum(2022,9,27,0,0,0):1/24:datenum(2022,10,2,0,0,0); %origin_date should Cover the Initial time
buffer_day = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------NO NEED TO CHANGE BELOW--------------------------
year = datevec(merge_date(1));
year = year(1);
outname = ['merged_forcing_',num2str(year),'.mat'];
for var_n = 1:length(merge_var)
    
    load(strcat(cfsr_dir,'\atm_forcing_',num2str(year),'.mat'),merge_var{var_n});
    load(strcat(cfsr_dir,'\atm_forcing_',num2str(year),'.mat'),'out_date');
    eval(strcat('cfsr_tmp = ',merge_var{var_n},';'));
    cfsr_date = out_date;

    variableInfo = who('-file', strcat('./wrf_forcing_',num2str(year),'.mat'));

    if(ismember(merge_var{var_n}, variableInfo))
        load(strcat('./wrf_forcing_',num2str(year),'.mat'),merge_var{var_n});
        load(strcat('./wrf_forcing_',num2str(year),'.mat'),'out_date');
        eval(strcat('wrf_tmp = ',merge_var{var_n},';'));
        wrf_date = out_date;
    else
        wrf_date= out_date;
        eval(strcat('wrf_tmp = ',merge_var{var_n},';'));
    end

    tmp = find(cfsr_date>=merge_date(1));
    merge_cfsr_start = tmp(1);
    tmp = find(cfsr_date<=merge_date(end));
    merge_cfsr_end = tmp(end);

    out_date = [cfsr_date(1:merge_cfsr_start-1), merge_date, cfsr_date(merge_cfsr_end+1:end)];

    buffer_pos_start = find(merge_date<=merge_date(1)+buffer_day);
    buffer_pos_end = find(merge_date>merge_date(end)-buffer_day);
    tmp = ones(1,length(merge_date));
    tmp(buffer_pos_start) = linspace(0,1,length(buffer_pos_start));
    tmp(buffer_pos_end) = linspace(1,0,length(buffer_pos_end));
    merge_coef = tmp;

    [r,c,~] = size(cfsr_tmp);

    for i = 1:r
        for j = 1:c
            tmp1 = interp1(cfsr_date,squeeze(cfsr_tmp(i,j,:)),merge_date,'linear').*(1-merge_coef);
            tmp2 = interp1(wrf_date,squeeze(wrf_tmp(i,j,:)),merge_date,'linear').*merge_coef;
            tmp0(i,j,:) = tmp1+tmp2;
        end
    end

    tmp = zeros(r,c,length(out_date));
    tmp(:,:,1:merge_cfsr_start-1) = cfsr_tmp(:,:,1:merge_cfsr_start-1);
    tmp(:,:,[1:size(tmp0,3)]+merge_cfsr_start-1) = tmp0;
    tmp(:,:,merge_cfsr_start+size(tmp0,3):end) = cfsr_tmp(:,:,merge_cfsr_end+1:end);

    eval(strcat(merge_var{var_n},' = tmp;'));
end

%Saving
if(~exist(outname,'file'))
    save(outname,"out_date",'-v7.3')
end

for var_n = 1:length(merge_var)
    eval(strcat('save(outname,','merge_var{var_n},''-append'');'));
end