function [ dat_final ] = read_woa_ini( fdir, lon, lat, sc, h, time, vname)
%   Detailed explanation goes here
    dis_c = 1;
    fname = dir([fdir '*.nc']);
    for i = 1:size(fname,1)
        filename{i,1} = strcat(fdir,fname(i).name);
    end
    [~,c] = size(filename);

    fnum = 14;
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
        case  {1 2 3}
            season = 1;
        case {4 5 6}
            season = 2;
        case {7 8 9}
            season = 3;
        case {10 11 12}
            season = 4;
    end

%     fnum = mon-1;
%     if(fnum<1)
%         fnum = mon+12;
%     end
%     fnum = fnum+1;

    fnum = 14;
    
    dat = ncread(filename{fnum},vname);
    dep = double(ncread(filename{fnum},'depth'));
    for j=1:k
        dat_pick(j,:) = squeeze(dat(i_pick(j),j_pick(j),:));
    end
    [~,~,z] = size(dat);

    for j=1:z
        tmp2 = dat_pick(:,j);
        tmp3 = griddata(lon_pick,lat_pick,tmp2,lon,lat); 
        try
            tmp3(isnan(tmp3)) = griddata(lon_pick(~isnan(tmp2)),lat_pick(~isnan(tmp2)),...
                tmp2(~isnan(tmp2)),lon(isnan(tmp3)),lat(isnan(tmp3)),'nearest');
            %dat_pick2(:,:,j) = tmp3.*mask;
            dat_pick2(:,:,j) = tmp3;
        catch 
            dat_pick2(:,:,j) = dat_pick2(:,:,j-1);
        end
            
    end
    

    [r,c,s] = size(sc);
    dat_final = zeros(r,c,s);
    for i=1:r
        for j=1:c
                raw_dep = dep;
                raw_dat = squeeze(dat_pick2(i,j,:));
                out_dep = squeeze(sc(i,j,:)).*h(i,j);
                tmp = interp1(raw_dep(~isnan(raw_dat)),raw_dat(~isnan(raw_dat)),out_dep,'spline');
                dat_final(i,j,:) = tmp;
        end
    end
    

        
end

