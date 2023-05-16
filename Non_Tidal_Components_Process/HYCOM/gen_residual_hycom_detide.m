% By Yuren Chen, 2021-11-26, Version 2.1
%--------------------------------------------------------------------------
% Program generate non-tidal boudaries from  output files of 
% Hybrid Coordinate Ocean Model (HYCOM).
% Due to the mismatch of land masks between hycom
% dataset and ROMS model grid, extrapolation are performed here. 
%--------------------------------------------------------------------------
clear all; close all;
for year = [2021]

% toolbox path
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\t_tide');
addpath(path,'../../../COAWST/Tools/mfiles/rutgers/utility');
addpath(path,'../../../COAWST/Tools/mfiles/roms_clm');

% Output time
date_out3d = datenum(year,1,1,0,0,0):6/24:datenum(year,12,31,24,0,0);
date_out2d = datenum(year,1,1,0,0,0):1/24:datenum(year,12,31,24,0,0);

% model grid
fn = '../../Model_grid/ROMS_WFS_10river_grid_v11.nc';

% params.
n_hycomlayer = 40;
%year = datevec(date_out3d(1));
%year = year(1);
N= length(ncread(fn,'Cs_r'));
n_layer = N;
Vtransform = ncread(fn,'Vtransform');                         
Vstretching = ncread(fn,'Vstretching');
THETA_S = ncread(fn,'theta_s');                     
THETA_B = ncread(fn,'theta_b');                     
TCLINE = ncread(fn,'Tcline');
hc = ncread(fn,'hc');
interp_scheme = 1; %1 for nearest, more in the future
dis_c = 0.1;
detide = 1;

%list for included hycom files
hycom_dir = ['../../HYCOM/',num2str(year),'/'];
fname = dir([hycom_dir 'HYCOM_*.nc']);
for i = 1:size(fname,1)
    filename{i,1} = strcat(hycom_dir,fname(i).name);
end
[n_days,c] = size(filename);

%--------------------- NO CHANGE BELOW-------------------------------------
%--------------------- Read from model grid--------------------------------
lon = ncread(fn,'lon_rho');
lat = ncread(fn,'lat_rho');
lon_u = ncread(fn,'lon_u');
lon_v = ncread(fn,'lon_v');
lat_u = ncread(fn,'lat_u');
lat_v = ncread(fn,'lat_v');
mask = ncread(fn,'mask_rho');
angle = ncread(fn,'angle');
h = ncread(fn,'h');

[r,c] = size(mask);

s_lon = lon(:,1:2);
s_lat = lat(:,1:2);
s_h = h(:,1:2);
s_mask = mask(:,1:2);
s_angle = angle(:,1:2);

n_lon = lon(:,end-1:end);
n_lat = lat(:,end-1:end);
n_h = h(:,end-1:end);
n_mask = mask(:,end-1:end);
n_angle = angle(:,end-1:end);

w_lon = lon(1:2,:);
w_lat = lat(1:2,:);
w_h = h(1:2,:);
w_mask = mask(1:2,:);
w_angle = angle(1:2,:);


e_lon = lon(end-1:end,:);
e_lat = lat(end-1:end,:);
e_h = h(end-1:end,:);
e_mask = mask(end-1:end,:);
e_angle = angle(end-1:end,:);

sb_num = length(s_lon(:,1));
nb_num = length(n_lon(:,1));
wb_num = length(w_lon(1,:));
eb_num = length(e_lon(1,:));

sb_range = [1:sb_num];
nb_range = [1:nb_num]+sb_range(end);
eb_range = [1:eb_num]+nb_range(end);
wb_range = [1:wb_num]+eb_range(end);

sb_range2 = [1:sb_num]+wb_range(end);
nb_range2 = [1:nb_num]+sb_range2(end);
eb_range2 = [1:eb_num]+nb_range2(end);
wb_range2 = [1:wb_num]+eb_range2(end);


% Calculate Sigma for each layer
sc = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,ones(size(h)),0,0);
sw = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,5,ones(size(h)),0,0);
sc = sc.*-1;
sw = sw.*-1;
dw = sw(:,:,1:end-1)-sw(:,:,2:end);
s_sc = squeeze(sc(:,1,:));
n_sc = squeeze(sc(:,end,:));
e_sc = squeeze(sc(end,:,:));
w_sc = squeeze(sc(1,:,:));

s_dw = squeeze(dw(:,1,:));
n_dw = squeeze(dw(:,end,:));
e_dw = squeeze(dw(end,:,:));
w_dw = squeeze(dw(1,:,:));

% all boundary points
b_lon = [s_lon(:,1)',n_lon(:,1)',e_lon(1,:),w_lon(1,:),s_lon(:,2)',n_lon(:,2)',e_lon(2,:),w_lon(2,:)];
b_lat = [s_lat(:,1)',n_lat(:,1)',e_lat(1,:),w_lat(1,:),s_lat(:,2)',n_lat(:,2)',e_lat(2,:),w_lat(2,:)];
b_sc = [s_sc;n_sc;e_sc;w_sc;s_sc;n_sc;e_sc;w_sc];
b_h = [s_h(:,1)',n_h(:,1)',e_h(1,:),w_h(1,:),s_h(:,2)',n_h(:,2)',e_h(2,:),w_h(2,:)];
b_mask = [s_mask(:,1)',n_mask(:,1)',e_mask(1,:),w_mask(1,:),s_mask(:,2)',n_mask(:,2)',e_mask(2,:),w_mask(2,:)];
b_angle = [s_angle(:,1)',n_angle(:,1)',e_angle(1,:),w_angle(1,:),s_angle(:,2)',n_angle(:,2)',e_angle(2,:),w_angle(2,:)];
b_num = length(b_lon);

day = 1;
depth = double(ncread(filename{day},'depth'));
lat_hycom = double(ncread(filename{day},'lat'));
lon_hycom = double(ncread(filename{day},'lon'));

% select data
[r_hycom,c_hycom] = size(lat_hycom);
k=0;
for i=1:r_hycom
    for j=1:c_hycom
        %dis = min(max([abs(lon_hycom(i,j)-b_lon);abs(lat_hycom(i,j)-b_lat)],[],1));
        dis = distance(lat_hycom(i,j),lon_hycom(i,j),b_lat,b_lon);
        if(min(abs(dis))<dis_c)
            k = k+1;
            i_pick(k) = i;
            j_pick(k) = j;
            lon_pick(k) = lon_hycom(i,j);
            lat_pick(k) = lat_hycom(i,j);
        end
    end
end

x0 = repmat(lon_pick',1,n_hycomlayer);
y0 = repmat(lat_pick',1,n_hycomlayer);
z0 = repmat(reshape(depth,1,n_hycomlayer),k,1);
x = repmat(b_lon',1,n_layer);
y = repmat(b_lat',1,n_layer);
z = squeeze(-1*set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,b_h',0,0));
figure;
scatter(b_lon,b_lat,40,'k','filled'); hold on;
scatter(lon_pick,lat_pick,40,'r','filled');hold on;

%--------------------- Read from model grid--------------------------------

data_time = [];
for day=1:n_days
    day
    tmp = double(ncread(filename{day},'time'))/24+datenum(2000,01,01);
    nrec = length(tmp);
    nrec_total = length(data_time);
    data_time = [data_time;tmp];   
    temp = double(ncread(filename{day},'temperature'));
    u = double(ncread(filename{day},'u'));
    v = double(ncread(filename{day},'v'));
    s = double(ncread(filename{day},'salinity'));
    el = double(ncread(filename{day},'surf_el'));
    
    clear el_pick s_pick temp_pick u_pick v_pick
    
    for i = 1:k
        u_pick(i,:,:) = u(i_pick(i),j_pick(i),:,:);
        v_pick(i,:,:) = v(i_pick(i),j_pick(i),:,:);
        el_pick(i,:) = el(i_pick(i),j_pick(i),:);
        temp_pick(i,:,:) = temp(i_pick(i),j_pick(i),:,:);
        s_pick(i,:,:) = s(i_pick(i),j_pick(i),:,:);
    end   
    
    for i=1:nrec
        el_t = el_pick(:,i);
 %       b_el(:,i+nrec_total) = barnes(lon_pick(~isnan(el_t)),lat_pick(~isnan(el_t)),el_t(~isnan(el_t)),b_lon,b_lat,1,1,3);
        b_el(:,i+nrec_total) = griddata(lon_pick(~isnan(el_t)),lat_pick(~isnan(el_t)),el_t(~isnan(el_t)),b_lon,b_lat,'linear');
        
        temp_t = temp_pick(:,:,i);
        s_t = s_pick(:,:,i);
        u_t = u_pick(:,:,i);
        v_t = v_pick(:,:,i);
    
        tmp = griddata(x0(~isnan(temp_t)),y0(~isnan(temp_t)),z0(~isnan(temp_t))...
            ,temp_t(~isnan(temp_t)),x,y,z,'linear');
        tmp(isnan(tmp)) = griddata(x0(~isnan(temp_t)),y0(~isnan(temp_t)),z0(~isnan(temp_t))...
            ,temp_t(~isnan(temp_t)),x(isnan(tmp)),y(isnan(tmp)),z(isnan(tmp)),'nearest');
 
        b_temp(:,:,i+nrec_total) = tmp;

            

        tmp = griddata(x0(~isnan(s_t)),y0(~isnan(s_t)),z0(~isnan(s_t))...
            ,s_t(~isnan(s_t)),x,y,z,'linear');
        tmp(isnan(tmp)) = griddata(x0(~isnan(s_t)),y0(~isnan(s_t)),z0(~isnan(s_t))...
            ,s_t(~isnan(s_t)),x(isnan(tmp)),y(isnan(tmp)),z(isnan(tmp)),'nearest');
 
        b_s(:,:,i+nrec_total) = tmp;

        b_u(:,:,i+nrec_total) = griddata(x0(~isnan(u_t)),y0(~isnan(u_t)),z0(~isnan(u_t)),u_t(~isnan(u_t)),x,y,z,'linear');
        b_v(:,:,i+nrec_total) = griddata(x0(~isnan(v_t)),y0(~isnan(v_t)),z0(~isnan(v_t)),v_t(~isnan(v_t)),x,y,z,'linear');
    end
end

if(detide)   
    b_u_tmp = tv3dinterpt(b_u,data_time,date_out2d,'linear',0);
    b_v_tmp = tv3dinterpt(b_v,data_time,date_out2d,'linear',0);
    b_el_tmp = tv2dinterpt(b_el,data_time,date_out2d,'linear',0);

    [a,b,~] = size(b_u_tmp);

    b_el2 = zeros(size(b_el));
    b_u2 = zeros(size(b_u));
    b_v2 = zeros(size(b_v));

    for i=1:a 
        if(b_mask(i)==1)
            if(sum(~isnan(b_el_tmp(i,:)))>0)
                [~,xout] = t_tide(b_el_tmp(i,:),'output','none');
                xout2 = interp1(date_out2d,xout,data_time);
                b_el2(i,:) = b_el(i,:)-xout2';
            end
            for j=1:b
                if(sum(~isnan(squeeze(b_u_tmp(i,j,:))))>0)
                    [~,xout] = t_tide(squeeze(b_u_tmp(i,j,:)),'output','none');
                    xout2 = interp1(date_out2d,xout,data_time);
                    b_u2(i,j,:) = squeeze(b_u(i,j,:))-xout2;

                    [~,xout] = t_tide(squeeze(b_v_tmp(i,j,:)),'output','none');
                    xout2 = interp1(date_out2d,xout,data_time);
                    b_v2(i,j,:) = squeeze(b_v(i,j,:))-xout2;
                end
            end
        end
    end
end

nrec_total = nrec_total+nrec;
b_s2 = b_s;
b_temp2 = b_temp;
b_el2(isnan(b_el2)) = 0;
b_u2(isnan(b_u2)) = 0;
b_v2(isnan(b_v2)) = 0;

[r,c] = size(lon);

%------------------interpolate to out date---------------------------------
b_temp3 = tv3dinterpt(b_temp2,data_time,date_out3d,'linear',1);
b_s3 = tv3dinterpt(b_s2,data_time,date_out3d,'linear',1);
b_u3 = tv3dinterpt(b_u2,data_time,date_out3d,'linear',1);
b_v3 = tv3dinterpt(b_v2,data_time,date_out3d,'linear',1);
b_el3 = tv2dinterpt(b_el2,data_time,date_out2d,'linear',1);
%------------------interpolate to out date---------------------------------       

%--------------------------rotate u,v to x,e-------------------------------


s_el3(:,1,:) = b_el3(sb_range,:);
n_el3(:,1,:) = b_el3(nb_range,:);
w_el3(:,1,:) = b_el3(wb_range,:);
e_el3(:,1,:) = b_el3(eb_range,:);

s_el3(:,2,:) = b_el3(sb_range2,:);
n_el3(:,2,:) = b_el3(nb_range2,:);
w_el3(:,2,:) = b_el3(wb_range2,:);
e_el3(:,2,:) = b_el3(eb_range2,:);

s_el3(isnan(s_el3)) = 0;
n_el3(isnan(n_el3)) = 0;
w_el3(isnan(w_el3)) = 0;
e_el3(isnan(e_el3)) = 0;

s_s3(:,:,:) = b_s3(sb_range,:,:);
n_s3(:,:,:) = b_s3(nb_range,:,:);
w_s3(:,:,:) = b_s3(wb_range,:,:);
e_s3(:,:,:) = b_s3(eb_range2,:,:);

s_temp3(:,:,:) = b_temp3(sb_range,:,:);
n_temp3(:,:,:) = b_temp3(nb_range,:,:);
w_temp3(:,:,:) = b_temp3(wb_range,:,:);
e_temp3(:,:,:) = b_temp3(eb_range2,:,:);

s_u3(:,:,1,:) = b_u3(sb_range,:,:);
n_u3(:,:,1,:) = b_u3(nb_range,:,:);
w_u3(:,:,1,:) = b_u3(wb_range,:,:);
e_u3(:,:,1,:) = b_u3(eb_range,:,:);

s_u3(:,:,2,:) = b_u3(sb_range2,:,:);
n_u3(:,:,2,:) = b_u3(nb_range2,:,:);
w_u3(:,:,2,:) = b_u3(wb_range2,:,:);
e_u3(:,:,2,:) = b_u3(eb_range2,:,:);

s_v3(:,:,1,:) = b_v3(sb_range,:,:);
n_v3(:,:,1,:) = b_v3(nb_range,:,:);
w_v3(:,:,1,:) = b_v3(wb_range,:,:);
e_v3(:,:,1,:) = b_v3(eb_range,:,:);

s_v3(:,:,2,:) = b_v3(sb_range2,:,:);
n_v3(:,:,2,:) = b_v3(nb_range2,:,:);
w_v3(:,:,2,:) = b_v3(wb_range2,:,:);
e_v3(:,:,2,:) = b_v3(eb_range2,:,:);

clear b_temp2 b_temp3 b_s2 b_s3 b_v2 b_v3 b_u2 b_u3

%rotate u,v to x,e

[~,~,~,b] = size(s_u3);


s_lon_u = rho2u_2d_mw(s_lon(:,1));
s_lon_v = rho2v_2d_mw(s_lon);
n_lon_u = rho2u_2d_mw(n_lon(:,1));
n_lon_v = rho2v_2d_mw(n_lon);
w_lon_u = rho2u_2d_mw(w_lon);
w_lon_v = rho2v_2d_mw(w_lon(1,:));
e_lon_u = rho2u_2d_mw(e_lon);
e_lon_v = rho2v_2d_mw(e_lon(2,:));

s_lat_u = rho2u_2d_mw(s_lat(:,1));
s_lat_v = rho2v_2d_mw(s_lat);
n_lat_u = rho2u_2d_mw(n_lat(:,1));
n_lat_v = rho2v_2d_mw(n_lat);
w_lat_u = rho2u_2d_mw(w_lat);
w_lat_v = rho2v_2d_mw(w_lat(1,:));
e_lat_u = rho2u_2d_mw(e_lat);
e_lat_v = rho2v_2d_mw(e_lat(2,:));

h_u = rho2u_2d_mw(h);
h_v = rho2v_2d_mw(h);

for i= 1:b
    for n_l= 1:n_layer
        s_u4(:,n_l,i) = rho2u_2d_mw(squeeze(s_u3(:,n_l,1,i)));
        n_u4(:,n_l,i) = rho2u_2d_mw(squeeze(n_u3(:,n_l,1,i)));
        w_u4(:,n_l,i) = rho2u_2d_mw([squeeze(w_u3(:,n_l,:,i))]');
        e_u4(:,n_l,i) = rho2u_2d_mw([squeeze(e_u3(:,n_l,:,i))]');

        s_v4(:,n_l,i) = rho2v_2d_mw(squeeze(s_v3(:,n_l,:,i)));
        n_v4(:,n_l,i) = rho2v_2d_mw(squeeze(n_v3(:,n_l,:,i)));
        w_v4(:,n_l,i) = rho2v_2d_mw([squeeze(w_v3(:,n_l,1,i))]');
        e_v4(:,n_l,i) = rho2v_2d_mw([squeeze(e_v3(:,n_l,2,i))]');

        s_uv3(:,n_l,i) = rho2u_2d_mw(squeeze(s_v3(:,n_l,1,i)));
        n_uv3(:,n_l,i) = rho2u_2d_mw(squeeze(n_v3(:,n_l,1,i)));
        w_uv3(:,n_l,i) = rho2u_2d_mw([squeeze(w_v3(:,n_l,:,i))]');
        e_uv3(:,n_l,i) = rho2u_2d_mw([squeeze(e_v3(:,n_l,:,i))]');

        s_vu3(:,n_l,i) = rho2v_2d_mw(squeeze(s_u3(:,n_l,:,i)));
        n_vu3(:,n_l,i) = rho2v_2d_mw(squeeze(n_u3(:,n_l,:,i)));
        w_vu3(:,n_l,i) = rho2v_2d_mw([squeeze(w_u3(:,n_l,1,i))]');
        e_vu3(:,n_l,i) = rho2v_2d_mw([squeeze(e_u3(:,n_l,2,i))]');  
    end    
end

s_angle_u = rho2u_2d_mw(s_angle(:,1));
n_angle_u = rho2u_2d_mw(n_angle(:,1));
w_angle_u = rho2u_2d_mw(w_angle);
e_angle_u = rho2u_2d_mw(e_angle);

s_angle_v = rho2v_2d_mw(s_angle);
n_angle_v = rho2v_2d_mw(n_angle);
w_angle_v = rho2v_2d_mw(w_angle(1,:));
e_angle_v = rho2v_2d_mw(e_angle(2,:));

s_mask_u = rho2u_2d_mw(s_mask(:,1));
n_mask_u = rho2u_2d_mw(n_mask(:,1));
w_mask_u = rho2u_2d_mw(w_mask);
e_mask_u = rho2u_2d_mw(e_mask);

s_mask_v = rho2v_2d_mw(s_mask);
n_mask_v = rho2v_2d_mw(n_mask);
w_mask_v = rho2v_2d_mw(w_mask(1,:));
e_mask_v = rho2v_2d_mw(e_mask(2,:));

s_el4(:,:) = squeeze(s_el3(:,1,:));
n_el4(:,:) = squeeze(n_el3(:,1,:));
w_el4(:,:) = squeeze(w_el3(:,1,:));
e_el4(:,:) = squeeze(e_el3(:,1,:));

for t=1:b
    for n_l = 1:n_layer
        s_u5(:,n_l,t) = cos(s_angle_u).*s_u4(:,n_l,t)+sin(s_angle_u).*s_uv3(:,n_l,t);
        s_v5(:,n_l,t) = -sin(s_angle_v).*s_vu3(:,n_l,t)+cos(s_angle_v).*s_v4(:,n_l,t);
        n_u5(:,n_l,t) = cos(n_angle_u).*n_u4(:,n_l,t)+sin(n_angle_u).*n_uv3(:,n_l,t);
        n_v5(:,n_l,t) = -sin(n_angle_v).*n_vu3(:,n_l,t)+cos(n_angle_v).*n_v4(:,n_l,t);
        e_u5(:,n_l,t) = cos(e_angle_u').*e_u4(:,n_l,t)+sin(e_angle_u').*e_uv3(:,n_l,t);
        e_v5(:,n_l,t) = -sin(e_angle_v').*e_vu3(:,n_l,t)+cos(e_angle_v').*e_v4(:,n_l,t);
        w_u5(:,n_l,t) = cos(w_angle_u').*w_u4(:,n_l,t)+sin(w_angle_u').*w_uv3(:,n_l,t);
        w_v5(:,n_l,t) = -sin(w_angle_v').*w_vu3(:,n_l,t)+cos(w_angle_v').*w_v4(:,n_l,t);

        s_u5(s_mask_u==0,t) = 0;
        s_v5(s_mask_v==0,t) = 0;

        n_u5(n_mask_u==0,t) = 0;
        n_v5(n_mask_v==0,t) = 0;

        e_u5(e_mask_u==0,t) = 0;
        e_v5(e_mask_v==0,t) = 0;

        w_u5(w_mask_u==0,t) = 0;
        w_v5(w_mask_v==0,t) = 0;
    end
end
%--------------------------rotate u,v to x,e-------------------------------

%---------------------------------2D  VELOCITY-----------------------------
e_dw_u = e_dw(:,:);
e_dw_v = e_dw(1:end-1,:);
w_dw_u = w_dw(:,:);
w_dw_v = w_dw(1:end-1,:);

s_dw_u = s_dw(1:end-1,:);
s_dw_v = s_dw(:,:);
n_dw_u = n_dw(1:end-1,:);
n_dw_v = n_dw(:,:);

for t = 1:length(date_out3d)
    e_u2d(:,t) = sum(e_u5(:,:,t).*e_dw_u,2);
    s_u2d(:,t) = sum(s_u5(:,:,t).*s_dw_u,2);
    w_u2d(:,t) = sum(w_u5(:,:,t).*w_dw_u,2);
    n_u2d(:,t) = sum(n_u5(:,:,t).*n_dw_u,2);
    
    e_v2d(:,t) = sum(e_v5(:,:,t).*e_dw_v,2);
    s_v2d(:,t) = sum(s_v5(:,:,t).*s_dw_v,2);
    w_v2d(:,t) = sum(w_v5(:,:,t).*w_dw_v,2);
    n_v2d(:,t) = sum(n_v5(:,:,t).*n_dw_v,2);
end

for i = 1:size(e_u2d,1)
    e_u2d_out(i,:) = interp1(date_out3d,e_u2d(i,:),date_out2d);
    w_u2d_out(i,:) = interp1(date_out3d,w_u2d(i,:),date_out2d);
end

for i = 1:size(s_u2d,1)
    s_u2d_out(i,:) = interp1(date_out3d,s_u2d(i,:),date_out2d);
    n_u2d_out(i,:) = interp1(date_out3d,n_u2d(i,:),date_out2d);
end

for i = 1:size(e_v2d,1)
    e_v2d_out(i,:) = interp1(date_out3d,e_v2d(i,:),date_out2d);
    w_v2d_out(i,:) = interp1(date_out3d,w_v2d(i,:),date_out2d);
end

for i = 1:size(s_v2d,1)
    s_v2d_out(i,:) = interp1(date_out3d,s_v2d(i,:),date_out2d);
    n_v2d_out(i,:) = interp1(date_out3d,n_v2d(i,:),date_out2d);
end
%Output

s_el_out = s_el4;
s_s_out = s_s3;
s_temp_out = s_temp3;
s_u_out = s_u5;
s_v_out = s_v5;

n_el_out = n_el4;
n_s_out = n_s3;
n_temp_out = n_temp3;
n_u_out = n_u5;
n_v_out = n_v5;

w_el_out = w_el4;
w_s_out = w_s3;
w_temp_out = w_temp3;
w_u_out = w_u5;
w_v_out = w_v5;

e_el_out = e_el4;
e_s_out = e_s3;
e_temp_out = e_temp3;
e_u_out = e_u5;
e_v_out = e_v5;

save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'date_out2d','date_out3d','-v7.3');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'s_el_out','-append');           
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'s_s_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'s_temp_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'s_u_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'s_v_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'s_u2d_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'s_v2d_out','-append');

save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'n_el_out','-append');           
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'n_s_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'n_temp_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'n_u_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'n_v_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'n_u2d_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'n_v2d_out','-append');

save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'w_el_out','-append');           
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'w_s_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'w_temp_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'w_u_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'w_v_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'w_u2d_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'w_v2d_out','-append');

save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'e_el_out','-append');           
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'e_s_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'e_temp_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'e_u_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'e_v_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'e_u2d_out','-append');
save(strcat('non_tidal_bnd_',num2str(year),'.mat'),'e_v2d_out','-append');

%Figures for checking
figure;
quiver(s_lon_u,s_lat_u,s_u4(:,1,1),s_v4(1:end-1,1,1),0,'r');
hold on;
quiver(s_lon_u,s_lat_u,s_u_out(:,1,1),s_v_out(1:end-1,1,1),0);
title('Soutern Boundary Velocity');
legend('Before rotation','After rotation');

figure;
quiver(n_lon_u,n_lat_u,n_u4(:,1,1),n_v4(1:end-1,1,1),0,'r');
hold on;
quiver(n_lon_u,n_lat_u,n_u_out(:,1,1),n_v_out(1:end-1,1,1),0);
title('Northern Boundary Velocity');
legend('Before rotation','After rotation');

figure;
quiver(w_lon_v',w_lat_v',w_u4(1:end-1,1,1),w_v4(:,1,1),0,'k');
hold on;
quiver(w_lon_v',w_lat_v',w_u_out(1:end-1,1,1),w_v_out(:,1,1),0,'g');
title('Western Boundary Velocity');
legend('Before rotation','After rotation');

figure;
quiver(e_lon_v',e_lat_v',e_u4(1:end-1,1,1),e_v4(:,1,1),0,'r');
hold on;
quiver(e_lon_v',e_lat_v',e_u_out(1:end-1,1,1),e_v_out(:,1,1),0);
title('Eastern Boundary Velocity');
legend('Before rotation','After rotation');

figure;
subplot(1,4,1);
contourf(w_temp3(:,:,1));
title('Western Boundary Temp.');
subplot(1,4,2);
contourf(e_temp3(:,:,1));
title('Eastern Boundary Temp.');
subplot(1,4,3);
contourf(s_temp3(:,:,1));
title('Soutern Boundary Temp.');
subplot(1,4,4);
contourf(n_temp3(:,:,1));
title('Northern Boundary Temp.');

figure;
subplot(1,4,1);
contourf(w_s3(:,:,20));
title('Western Boundary Sal.');
subplot(1,4,2);
contourf(e_s3(:,:,20));
title('Eastern Boundary Sal.');
subplot(1,4,3);
contourf(s_s3(:,:,20));
title('Soutern Boundary Sal.');
subplot(1,4,4);
contourf(n_s3(:,:,20));
title('Northern Boundary Sal.');

clear all; close all;
end



