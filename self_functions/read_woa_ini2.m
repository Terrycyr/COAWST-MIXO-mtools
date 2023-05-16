function [ dat_final, dat_m, dep ] = read_woa_ini2( fdir, lon, lat, sc, h, time, vname)
%   Detailed explanation goes here
    dis_c = 1;
    %dep_range = [40 200];
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
        case {1 2 3}
            season = 1;
        case {4 5 6}
            season = 2;
        case {7 8 9}
            season = 3;
        case {10 11 12}
            season = 4;
    end

    fnum = 1;
    dep_mon = double(ncread(filename{fnum},'depth'));
    dat_mon = ncread(filename{fnum},vname);

    dat = dat_mon;
    dep = dep_mon;

    for j=1:k
        dat_pick(j,:) = squeeze(dat(i_pick(j),j_pick(j),:));
    end
    [~,~,z] = size(dat);
    
    [~,~,N] = size(sc);

    dat_m = mean(dat_pick,1,'omitnan');
    %dat_m = smooth(dat_m,0.2,'loess');
    dat_m(dat_m<0) = 0;
    dat_m(isnan(dat_m)) = 0;
    
    %pos1(1)=sum(dep<=dep_range(1));
    %pos1(2)=sum(dep<=dep_range(2));

    %scale = exp(linspace(-16,0,length(pos1(1):pos1(2))));
    %k=0;
    %for i=pos1(1):pos1(2)
    %    k=k+1;
    %    dat_m(i) = dat_m(pos1(1))*(1-scale(k))+dat_m(pos1(2))*scale(k);
    %end

    figure;
    plot(dep,dat_m);
    
    [r,c] = size(lon);
    for i=1:r
        for j=1:c
            for k=1:N
                dat_final(i,j,k) = interp1(dep, dat_m,h(i,j).*sc(i,j,k),'linear');
            end
        end
    end
end

