function [u_out, v_out, u4, v4] = get_OTPS_uv_result(grd,out_time,OTPS_result_u,OTPS_result_v)
%get_OTPS_uv_result Summary of this function goes here
%   Detailed explanation goes here

%--------------------- Read from model grid--------------------------------
    fid = fopen(OTPS_result_u);
    fid2 = fopen(OTPS_result_v);
    for i=1:6
       fgetl(fid);
       fgetl(fid2);
    end
   
    for i=1:length(grd.blon_u)
        for t=1:length(out_time)
            dat=fgetl(fid);
            if(strcmp(dat(22),'*'))
                u(i,t) = 0;
                uv(i,t) = 0;
                dep_u(i) = grd.bh_u(i);
            else
                u(i,t) = str2num(dat(66:76))/100;
                uv(i,t) = str2num(dat(77:86))/100;
                dep_u(i) = str2num(dat(87:end));
            end 
        end
    end
    
    for i=1:length(grd.blon_v)
        for t=1:length(out_time)
            dat=fgetl(fid2);
            if(strcmp(dat(22),'*'))
                v(i,t) = 0;
                vu(i,t) = 0;
                dep_v(i) = grd.bh_v(i);
            else
                v(i,t) = str2num(dat(77:86))/100;
                vu(i,t) = str2num(dat(66:76))/100;
                dep_v(i) = str2num(dat(87:end));
            end 
        end  
    end
   
    
    % Find mismatch wet points 
    for i=2:length(grd.blon_u)
        dis_u(i-1) = distance(grd.blat_u(1),grd.blon_u(1),grd.blat_u(i),grd.blon_u(i));
    end
    
    for i=2:length(grd.blon_v)
        dis_v(i-1) = distance(grd.blat_v(1),grd.blon_v(1),grd.blat_v(i),grd.blon_v(i));
    end
    
    dis_u2 = [0,dis_u];
    dis_v2 = [0,dis_v];


    mis_u = [sum([u>0],2)>0] + grd.bmask_u;
    mis_v = [sum([v>0],2)>0] + grd.bmask_v;

    mis_u(grd.bmask_u==0) =0;
    mis_v(grd.bmask_v==0) =0;


    % Use linear estrapolation to eliminate the mismatch
    for t=1:length(out_time)
        u(mis_u==1,t) = interp1(dis_u2(mis_u==2),u(mis_u==2,t),...
                dis_u2(mis_u==1),'linear','extrap');
        uv(mis_u==1,t) = interp1(dis_u2(mis_u==2),uv(mis_u==2,t),...
                dis_u2(mis_u==1),'linear','extrap');
        v(mis_v==1,t) = interp1(dis_v2(mis_v==2),v(mis_v==2,t),...
                dis_v2(mis_v==1),'linear','extrap');
        vu(mis_v==1,t) = interp1(dis_v2(mis_v==2),vu(mis_v==2,t),...
                dis_v2(mis_v==1),'linear','extrap');
    end

    %rotate u,v to x,e
    u4 = u;
    v4 = v;


    [a,b] = size(u4);

    for t=1:b   
        u(:,t) = cos(grd.bangle_u).*u4(:,t)+sin(grd.bangle_u).*uv(:,t);
        v(:,t) = -sin(grd.bangle_v).*vu(:,t)+cos(grd.bangle_v).*v4(:,t);

    % Flux Correction
        if(t==1)
            cor_u = dep_u'./grd.bh_u;
            cor_u(grd.bh_u<20) = 1;
            cor_v = dep_v'./grd.bh_v;
            cor_v(grd.bh_v<20) = 1;
        end
        u(:,t) = u(:,t).*cor_u;
        v(:,t) = v(:,t).*cor_v;
        
        u(grd.bmask_u==0,t) = 0;
        v(grd.bmask_v==0,t) = 0;
    end
    
    u_out = u;
    v_out = v;
    
    figure;
    plot(-1*dep_u,'r');hold on;
    plot(-1*grd.bh_u,'k');
    legend('OTPS','MODEL');
end

