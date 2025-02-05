function [ i,j ] = find_ij( flon, flat, lon, lat, mask)
%find_ij Summary of this function goes here
%   Detailed explanation goes here
    if(flon<=max(max(lon))&&flon>=min(min(lon))&&flat<=max(max(lat))&&flat>=min(min(lat)))
        dis = distance(flat,flon,lat,lon);
        dis(mask==0) = 1000000;
        [i,j] = find(dis==min(min(dis)));
    else
        'Out of area!'
    end
end

