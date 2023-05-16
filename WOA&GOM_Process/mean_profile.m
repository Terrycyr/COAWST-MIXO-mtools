    clear all; close all;
    addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');
    
    dis_c = 1;
    fname = dir([fdir '*.nc']);
    for i = 1:size(fname,1)
        filename{i,1} = strcat(fdir,fname(i).name);
    end
    [~,c] = size(filename);

    fnum = 1;
    depth = double(ncread(filename{fnum},'depth'));
    lat_dataset = double(ncread(filename{fnum},'lat'));
    lon_dataset = double(ncread(filename{fnum},'lon'));
    [lat_dataset,lon_dataset] = meshgrid(lat_dataset,lon_dataset);

    % select data
    [r_dataset,c_dataset] = size(lat_dataset);
    k=0;
    for i=1:r_dataset
        for j=1:c_dataset
            dis = min(max([abs(lon_dataset(i,j)-lon(:)');abs(lat_dataset(i,j)-lat(:)')],[],1));
            if(min(abs(dis))<dis_c)
                k = k+1;
                i_pick(k) = i;
                j_pick(k) = j;
                lon_pick(k) = lon_dataset(i,j);
                lat_pick(k) = lat_dataset(i,j);
            end
        end
    end

    figure;
    scatter(lon(:),lat(:),40,'k','filled'); hold on;
    scatter(lon_pick,lat_pick,40,'r','filled');hold on;

    tvec = datevec(time);
    mon = tvec(2);

    switch mon
        case {12 1 2}
            season = 4;
        case {3 4 5}
            season = 1;
        case {6 7 8}
            season = 2;
        case {9 10 11}
            season = 3;
    end

    dep_mon = double(ncread(filename{mon},'depth'));
    dat_mon = ncread(filename{mon},vname);
    dat = dat_mon;
    dep = dep_mon;

    for j=1:k
        dat_pick(j,:) = squeeze(dat(i_pick(j),j_pick(j),:));
    end
    [~,~,z] = size(dat);
    
    [~,~,N] = size(sc);

    dat_m = median(dat_pick,1,'omitnan');
    %dat_std = std(dat_pick,'omitnan');
    %dat_std(isnan(dat_std)) = 0;
    %dat_m = dat_m-dat_std;
    dat_m = smooth(dat_m,0.2,'loess');
    dat_m(dat_m<0) = 0;
    dat_m(isnan(dat_m)) = 0;