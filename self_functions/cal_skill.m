function [skill] = cal_skill(m,obs)
%cal_skill Summary of this function goes here
%   Detailed explanation goes here
    d1 = (m-obs).^2;
    d2 = (obs-mean(obs,'omitnan')).^2;
    d1(isnan(d1)) = 0;
    d2(isnan(d2)) = 0;
    skill = 1-(sum(d1)/sum(d2));
end

