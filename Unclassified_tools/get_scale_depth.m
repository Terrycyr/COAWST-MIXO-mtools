function [scale] = get_scale_depth(dep,bay_mask,all_mask,lb,snum)
%get_scale_distance Summary of this function goes here
%   Detailed explanation goes here

[r,c] = size(all_mask);
dep(all_mask==0)=NaN;
depmin = min(min(dep,[],'omitnan'));
scale = zeros(r,c);

for i=1:r
    for j=1:c
        if(bay_mask(i,j)==0&&all_mask(i,j)==1)
            scale(i,j) = snum^(dep(i,j)-depmin);
            if(scale(i,j)<lb)
                scale(i,j)=lb;
            end
        elseif(bay_mask(i,j)==1)
            scale(i,j) = 1;
        end
    end
end

end