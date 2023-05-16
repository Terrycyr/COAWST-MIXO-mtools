function [el_out] = get_OTPS_el_result(grd,out_time,OTPS_result)
%get_OTPS_result Summary of this function goes here
%   Detailed explanation goes here
%--------------------- Read from model grid--------------------------------
    fid = fopen(OTPS_result);
    
% Skip header lines
    for i=1:6
       fgetl(fid);
    end

    for i=1:length(grd.blon)
        for t=1:length(out_time)
            dat=fgetl(fid);
            if(strcmp(dat(22),'*'))
                el(i,t) = 0;
                dep(i) = grd.bh(i);
            else
                el(i,t) = str2num(dat(46:56));
                dep(i) = str2num(dat(57:end));
            end
        end   
    end

    for i=2:length(grd.blon)
        dis(i-1) = distance(grd.blat(1),grd.blon(1),grd.blat(i),grd.blon(i));
    end

    dis2 = [0,dis];
    mis = double([sum([el>0],2)>0]) + grd.bmask;
    mis(grd.bmask==0) =0;

    for t=1:length(out_time)
        el(mis==1,t) = interp1(dis2(mis==2),el(mis==2,t),...
                dis2(mis==1),'linear','extrap');
        el(:,t) = el(:,t).*grd.bmask;
    end

    el_out = el;
    
end

