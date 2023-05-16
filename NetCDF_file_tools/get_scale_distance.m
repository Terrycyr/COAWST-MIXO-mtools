function [scale] = get_scale_distance(lon,lat,bay_mask,all_mask,lb,snum)
%get_scale_distance Summary of this function goes here
%   Detailed explanation goes here

[r,c] = size(all_mask);
scale = zeros(r,c);

for i=1:r
    for j=1:c
        if(bay_mask(i,j)==0&&all_mask(i,j)==1)
            dis = distance(lat(i,j),lon(i,j),lat(bay_mask==1),lon(bay_mask==1));
            dis2 = min(dis);
            scale(i,j) = snum^dis2;
            if(scale(i,j)<lb)
                scale(i,j)=lb;
            end
        elseif(bay_mask(i,j)==1)
            scale(i,j) = 1;
        end
    end
end

end