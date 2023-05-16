function [ dat_out ] = read_woa( fdir, b_lon, b_lat, b_mask, b_sc, b_h, time, vname)
%   Detailed explanation goes here
    dis_c = 1;
    fname = dir([fdir '*.nc']);
    for i = 1:size(fname,1)
        filename{i,1} = strcat(fdir,fname(i).name);
    end
    [~,c] = size(filename);

    fnum = 2;
    depth = double(ncread(filename{fnum},'depth'));
    lat_dataset = double(ncread(filename{fnum},'lat'));
    lon_dataset = double(ncread(filename{fnum},'lon'));
    [lat_dataset,lon_dataset] = meshgrid(lat_dataset,lon_dataset);

    % select data
    [r_dataset,c_dataset] = size(lat_dataset);
    k=0;
    for i=1:r_dataset
        for j=1:c_dataset
            dis = min(max([abs(lon_dataset(i,j)-b_lon);abs(lat_dataset(i,j)-b_lat)],[],1));
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
    scatter(b_lon,b_lat,40,'k','filled'); hold on;
    scatter(lon_pick,lat_pick,40,'r','filled');hold on;

    i=0;
    for mon = 1:12
        i=i+1;
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
       
        
        fnum = season+12+1;
        dep_mon = double(ncread(filename{fnum},'depth'));
        dat_mon = ncread(filename{fnum},vname);

        dat = dat_mon;
        dep = dep_mon;

        for j=1:k
            dat_pick(j,:,i) = squeeze(dat(i_pick(j),j_pick(j),:));
        end
        [~,~,tmp] = size(dat);

        for j=1:tmp
            tmp2 = dat_pick(:,j,i);
            tmp3 = griddata(lon_pick,lat_pick,tmp2,b_lon,b_lat);
            try
                tmp3(isnan(tmp3)) = griddata(lon_pick(~isnan(tmp2)),lat_pick(~isnan(tmp2)),...
                    tmp2(~isnan(tmp2)),b_lon(isnan(tmp3)),b_lat(isnan(tmp3)),'nearest');
                dat_pick2(:,j) = tmp3.*b_mask;
            catch
                dat_pick2(:,j) = dat_pick2(:,j-1);
            end
        end

        for j=1:length(b_mask)
            dat_final(j,:,i) = interp1(dep,dat_pick2(j,:),b_sc(j,:).*b_h(j),'linear');
        end
    end
    
    i = 0;
    for t = time
        i = i+1;
        tvec = datevec(t);
        pos = tvec(2);
        dat_out(:,:,i) = dat_final(:,:,pos);
    end
        
end

