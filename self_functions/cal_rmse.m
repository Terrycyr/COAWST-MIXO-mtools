function [rmse] = cal_rmse(m,obs)
%cal_skill Summary of this function goes here
%   Detailed explanation goes here
    rmse = sqrt(mean((m-obs).^2,'omitnan'));   
end

