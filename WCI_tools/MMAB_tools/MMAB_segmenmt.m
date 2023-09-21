clear all
b_points = load('MMAB_boudary_point.dat');
% lat = load('/home/terry/SWAN/grid/grid/PRE_lat.dat');
% lon = load('/home/terry/SWAN/grid/grid/PRE_lon.dat');
% dep = load('/home/terry/SWAN/grid/PRE_depth.dat');
lat = load('../COAWST/pre_coawstswan_2017/PRE_2017_lat.dat');
lon = load('../COAWST/pre_coawstswan_2017/PRE_2017_lon.dat');
dep = load('../COAWST/pre_coawstswan_2017/PRE_2017_depth.dat');
mask = ncread('../COAWST/pre_coawstswan_2017/pre_2017_grd.nc','mask_rho');
mask = mask';
w_num = 4;
e_num = 4;
s_num = 9;
w_start = 1;
w_end = 76;
e_start = 1;
e_end = 95;
s_start = 1;
s_end = 627;
nx = 627;
ny = 546;
[r,c] = size(lat);
year = 2017
% w
for i=1:w_num
    d_old=1000;
    if(i==1) 
        for j=1:r
            d = distance(b_points(i,4),b_points(i,3),lat(j,1),lon(j,1));
            if(d<d_old)
                w_g(i) = j;
                d_old = d;
            end
        end
        w_seg_head(i) = w_start;
        w_d(i) = distance(lat(w_g(i),1),lon(w_g(i),1),lat(w_seg_head(i),1),lon(w_seg_head(i),1));
        if(w_num==1)
            w_seg_head(i+1) = w_end;
        end
    elseif(i<w_num)
        for j=1:r
            d = distance(b_points(i,4),b_points(i,3),lat(j,1),lon(j,1));
            if(d<d_old)
                w_g(i) = j;
                d_old = d;
            end
        end
        w_seg_head(i) = round((w_g(i)+w_g(i-1))/2);
        w_d(i) = distance(lat(w_g(i),1),lon(w_g(i),1),lat(w_seg_head(i),1),lon(w_seg_head(i),1));
    else
        for j=1:r
            d = distance(b_points(i,4),b_points(i,3),lat(j,1),lon(j,1));
            if(d<d_old)
                w_g(i) = j;
                d_old = d;
            end
        end
        w_seg_head(i) = round((w_g(i)+w_g(i-1))/2);
        w_d(i) = distance(lat(w_g(i),1),lon(w_g(i),1),lat(w_seg_head(i),1),lon(w_seg_head(i),1));
        w_seg_head(i+1) = w_end;
    end
end


%e
for i=1:e_num
    d_old=1000;
    if(i==1) 
        for j=1:r
            d = distance(b_points(i+w_num,4),b_points(i+w_num,3),lat(j,end),lon(j,end));
            if(d<d_old)
                e_g(i) = j;
                d_old = d;
            end
        end
        e_seg_head(i) = e_start;
        e_d(i) = distance(lat(e_g(i),end),lon(e_g(i),end),lat(e_seg_head(i),end),lon(e_seg_head(i),end));
        if(e_num==1)
            e_seg_head(i+1) = e_end;
        end
    elseif(i<e_num)
        for j=1:r
            d = distance(b_points(i+w_num,4),b_points(i+w_num,3),lat(j,end),lon(j,end));
            if(d<d_old)
                e_g(i) = j;
                d_old = d;
            end
        end
        e_seg_head(i) = round((e_g(i)+e_g(i-1))/2);
        e_d(i) = distance(lat(e_g(i),end),lon(e_g(i),end),lat(e_seg_head(i),end),lon(e_seg_head(i),end));
    else
        for j=1:r
            d = distance(b_points(i+w_num,4),b_points(i+w_num,3),lat(j,end),lon(j,end));
            if(d<d_old)
                e_g(i) = j;
                d_old = d;
            end
        end
        e_seg_head(i) = round((e_g(i)+e_g(i-1))/2);
        e_d(i) = distance(lat(e_g(i),end),lon(e_g(i),end),lat(e_seg_head(i),end),lon(e_seg_head(i),end));
        e_seg_head(i+1) = e_end;
    end
end
%s

for i=1:s_num
    d_old=1000;
    if(i==1) 
        for j=1:c
            d = distance(b_points(i+w_num+e_num,4),b_points(i+w_num+e_num,3),lat(1,j),lon(1,j));
            if(d<d_old)
                s_g(i) = j;
                d_old = d;
            end
        end
        s_seg_head(i) = s_start;
        s_d(i) = distance(lat(1,s_g(i)),lon(1,s_g(i)),lat(1,s_seg_head(i)),lon(1,s_seg_head(i)));
        if(s_num==1)
            s_seg_head(i+1) = s_end;
        end
    elseif(i<s_num)
        for j=1:c
            d = distance(b_points(i+w_num+e_num,4),b_points(i+w_num+e_num,3),lat(1,j),lon(1,j));
            if(d<d_old)
                s_g(i) = j;
                d_old = d;
            end
        end
        s_seg_head(i) = round((s_g(i)+s_g(i-1))/2);
        s_d(i) = distance(lat(1,s_g(i)),lon(1,s_g(i)),lat(1,s_seg_head(i)),lon(1,s_seg_head(i)));
    else
        for j=1:c
            d = distance(b_points(i+w_num+e_num,4),b_points(i+w_num+e_num,3),lat(1,j),lon(1,j));
            if(d<d_old)
                s_g(i) = j;
                d_old = d;
            end
        end
        s_seg_head(i) = round((s_g(i)+s_g(i-1))/2);
        s_d(i) = distance(lat(1,s_g(i)),lon(1,s_g(i)),lat(1,s_seg_head(i)),lon(1,s_seg_head(i)));
        s_seg_head(i+1) = s_end;
    end
end

fid = fopen('MMAB_segment.dat','wt+');
for i=1:w_num
    fprintf(fid,['BOUNdspec SEGment IJ 0 %d 0 %d VARiable FILE %8.4f ''MMAB_WAVE_',num2str(year),'_%d.bnd'' 1\n'],...
        w_seg_head(i)-1,w_seg_head(i+1)-1,w_d(i),i);
end

for i=1:e_num
    fprintf(fid,['BOUNdspec SEGment IJ ',num2str(nx-1),' %d ', num2str(nx-1),' %d VARiable FILE %8.4f ''MMAB_WAVE_',num2str(year),'_%d.bnd'' 1\n'],...
        e_seg_head(i)-1,e_seg_head(i+1)-1,e_d(i),i+w_num);
end

for i=1:s_num
    fprintf(fid,['BOUNdspec SEGment IJ %d 0 %d 0 VARiable FILE %8.4f ''MMAB_WAVE_',num2str(year),'_%d.bnd'' 1\n'],...
        s_seg_head(i)-1,s_seg_head(i+1)-1,s_d(i),i+w_num+e_num);
end

for i=1:w_num
    fprintf(fid,'POINTS ''w_%d''%8.3f%8.3f\n',i,lon(w_g(i),1),lat(w_g(i),1));
end

for i=1:e_num
    fprintf(fid,'POINTS ''e_%d''%8.3f%8.3f\n',i,lon(e_g(i),end),lat(e_g(i),end));
end

for i=1:s_num
    fprintf(fid,'POINTS ''s_%d''%8.3f%8.3f\n',i,lon(1,s_g(i)),lat(1,s_g(i)));
end