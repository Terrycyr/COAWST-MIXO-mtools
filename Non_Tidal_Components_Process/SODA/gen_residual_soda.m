% By Yuren Chen, 2020-07-15, Version 1.0
%--------------------------------------------------------------------------
% Program generate non-tidal boudaries from  output files of 
% Simple Ocean Data Assimilation (SODA).
% Due to the mismatch of land masks and bathymetry between SODA
% dataset and ROMS model grid, extrapolation and flux correction are
% performed here. 
%--------------------------------------------------------------------------

clear all; close all;

% toolbox path
addpath(path,'/home/terry/matlab_tool');
addpath(path,'/home/terry/COAWST/Tools/mfiles/rutgers/utility');
addpath(path,'/home/terry/COAWST/Tools/mfiles/roms_clm');

% Output time
date_out3d = datenum(2017,7,1,0,0,0):1:datenum(2017,9,17,24,0,0);
date_out2d = datenum(2017,7,1,0,0,0):1/24:datenum(2017,9,17,24,0,0);

% model grid
fn = '../../pre_new_grd.nc';

% params.
n_layer = 15;
n_sodalayer = 50;
n_sodaperd = 1;
year = 2017;
N= 15;
Vtransform = 2;                          
Vstretching = 4;
THETA_S = 3.0d0;                     
THETA_B = 3.0d0;                     
TCLINE = 0.0d0;
hc = TCLINE;
interp_scheme = 1; %1 for nearest, more in the future
dis_c = 2;

%list for included SODA files
filename = textread(['./list',num2str(year),'.dat'],'%s');
[n_days,c] = size(filename);

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
wb_num = length(w_lon(1,:));
eb_num = length(e_lon(1,:));

sb_range = [1:sb_num];
wb_range = [sb_num+1:sb_num+wb_num];
eb_range = [sb_num+wb_num+1:sb_num+wb_num+eb_num];

sb_range2 = [1:sb_num]+sb_num+wb_num+eb_num;
wb_range2 = [sb_num+1:sb_num+wb_num]+sb_num+wb_num+eb_num;
eb_range2 = [sb_num+wb_num+1:sb_num+wb_num+eb_num]+sb_num+wb_num+eb_num;


% Calculate Sigma for each layer
zlev = set_depth(Vtransform,Vstretching,THETA_S,THETA_B,hc,N,1,h,0,0);
sc = zlev./h;
sc = sc.*-1;
s_sc = squeeze(sc(:,1,:));
e_sc = squeeze(sc(end,:,:));
w_sc = squeeze(sc(1,:,:));

% all boundary points
b_lon = [s_lon(:,1)',e_lon(1,:),w_lon(1,:),s_lon(:,2)',e_lon(2,:),w_lon(2,:)];
b_lat = [s_lat(:,1)',e_lat(1,:),w_lat(1,:),s_lat(:,2)',e_lat(2,:),w_lat(2,:)];
b_sc = [s_sc;e_sc;w_sc;s_sc;e_sc;w_sc];
b_h = [s_h(:,1)',e_h(1,:),w_h(1,:),s_h(:,2)',e_h(2,:),w_h(2,:)];
b_mask = [s_mask(:,1)',e_mask(1,:),w_mask(1,:),s_mask(:,2)',e_mask(2,:),w_mask(2,:)];
b_angle = [s_angle(:,1)',e_angle(1,:),w_angle(1,:),s_angle(:,2)',e_angle(2,:),w_angle(2,:)];
b_num = 2*(r+2*c);

day = 1;
depth = double(ncread(filename{day},'st_ocean'));
lat_soda0 = double(ncread(filename{day},'yu_ocean'));
lon_soda0 = double(ncread(filename{day},'xu_ocean'));
lat_soda_t0 = double(ncread(filename{day},'yt_ocean'));
lon_soda_t0 = double(ncread(filename{day},'xt_ocean'));
[lat_soda,lon_soda] = meshgrid(lat_soda0,lon_soda0);
[lat_soda_t,lon_soda_t] = meshgrid(lat_soda_t0,lon_soda_t0);

% select data
[r_soda,c_soda] = size(lat_soda);
[r_soda_t,c_soda_t] = size(lat_soda_t);
k=0;
for i=1:r_soda
    for j=1:c_soda
        dis = distance(lat_soda(i,j),lon_soda(i,j),b_lat,b_lon);
        if(min(abs(dis))<dis_c)
            k = k+1;
            i_pick(k) = i;
            j_pick(k) = j;
            lon_pick(k) = lon_soda(i,j);
            lat_pick(k) = lat_soda(i,j);
        end
    end
end

k2=0;
for i=1:r_soda_t
    for j=1:c_soda_t
        dis = distance(lat_soda_t(i,j),lon_soda_t(i,j),b_lat,b_lon);
        if(min(abs(dis))<dis_c)
            k2 = k2+1;
            i_pick_t(k2) = i;
            j_pick_t(k2) = j;
            lon_pick_t(k2) = lon_soda_t(i,j);
            lat_pick_t(k2) = lat_soda_t(i,j);
        end
    end
end

figure;
scatter(b_lon,b_lat,40,'k','filled'); hold on;
scatter(lon_pick,lat_pick,40,'r','filled');hold on;
scatter(lon_pick_t,lat_pick_t,40,'g','filled');hold on;

%--------------------- Read from model grid--------------------------------

%-------------------------horizontal interplation--------------------------
% Barnes method is used in the horizontal interpolation and extrapolation
% of temp., sal. and zeta, while that of u and v is linear interoplation to
% avoid introducing extra flux. 

for day=1:n_days
    day
    data_day(day) = double(ncread(filename{day},'time'))+datenum(1980,01,01);
    
    temp = double(ncread(filename{day},'temp'));
    u = double(ncread(filename{day},'u'));
    v = double(ncread(filename{day},'v'));
    s = double(ncread(filename{day},'salt'));
    el = double(ncread(filename{day},'ssh'));
    
    clear el_pick s_pick temp_pick u_pick v_pick
    
    for i = 1:k
        u_pick(i,:) = u(i_pick(i),j_pick(i),:);
        v_pick(i,:) = v(i_pick(i),j_pick(i),:);
    end
    
    for i = 1:k2
        el_pick(i,:) = el(i_pick_t(i),j_pick_t(i),:);
        temp_pick(i,:) = temp(i_pick_t(i),j_pick_t(i),:);
        s_pick(i,:) = s(i_pick_t(i),j_pick_t(i),:);
    end
        
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
% Add in 0 m layer

        [a,b] = size(temp_pick);
        temp_t = [];
        temp_t(:,1) = temp_pick(:,1);
        temp_t(:,2:b+1) = temp_pick;
        temp_pick = temp_t;

        [a,b] = size(s_pick);
        s_t = [];
        s_t(:,1) = s_pick(:,1);
        s_t(:,2:b+1) = s_pick;
        s_pick = s_t;

        [a,b] = size(u_pick);
        u_t = [];
        u_t(:,1) = u_pick(:,1);
        u_t(:,2:b+1) = u_pick;
        u_pick = u_t;

        [a,b] = size(v_pick);
        v_t = [];
        v_t(:,1) = v_pick(:,1);
        v_t(:,2:b+1) = v_pick;
        v_pick = v_t;

    if(day==1)
        depth = [0;depth];
        n_sodalayer = n_sodalayer+1;
    end
    
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    for i=1:n_sodaperd
        b_el(:,i+(day-1)*n_sodaperd) = barnes(lon_pick_t(~isnan(el_pick)),lat_pick_t(~isnan(el_pick)),el_pick(~isnan(el_pick)),b_lon,b_lat,1,1,3);
        for j=1:n_sodalayer
            temp_t = temp_pick(:,j);
            s_t = s_pick(:,j);
            u_t = u_pick(:,j);
            v_t = v_pick(:,j);
            try
                b_temp(:,j,i+(day-1)*n_sodaperd) = barnes(lon_pick_t(~isnan(temp_t)),lat_pick_t(~isnan(temp_t)),temp_t(~isnan(temp_t)),...
                    b_lon(:),b_lat(:),1,1,3);
            catch
                b_temp(:,j,i+(day-1)*n_sodaperd) = 0;
            end
            
            try
                b_s(:,j,i+(day-1)*n_sodaperd) = barnes(lon_pick_t(~isnan(s_t)),lat_pick_t(~isnan(s_t)),s_t(~isnan(s_t)),...
                    b_lon(:),b_lat(:),1,1,3);
            catch
                b_s(:,j,i+(day-1)*n_sodaperd) = 0;
            end
            
            b_u(:,j,i+(day-1)*n_sodaperd) = griddata(lon_pick,lat_pick,u_t,b_lon(:),b_lat(:),'linear');
            b_v(:,j,i+(day-1)*n_sodaperd) = griddata(lon_pick,lat_pick,v_t,b_lon(:),b_lat(:),'linear');
        end
    end
end

b_u(isnan(b_u)) = 0;
b_v(isnan(b_v)) = 0;
%-------------------------horizontal interplation--------------------------



%-------------------------integrate to 2d----------------------------------
[r,c] = size(lon);

b_temp2 = zeros(r,n_layer,n_days*n_sodaperd);
b_s2 = zeros(r,n_layer,n_days*n_sodaperd);
b_u2 = zeros(r,n_days*n_sodaperd);
b_v2 = zeros(r,n_days*n_sodaperd);



for t = 1:n_days*n_sodaperd
    
%-------------------------interpolate to s-level---------------------------
% The method use for temp. and salinity is different from that for u and v.
% Vetically linear interpolation are applied to the former. 
% This can help to remain the vetical distribution of both active tracers. 
% And the latter is just adding up for integration, which mainly ensures the flux conservation.

    for i=1:b_num
        b_temp2(i,:,t) = interp1(depth,b_temp(i,:,t),b_h(i).*b_sc(i,:));
        
        b_s2(i,:,t) = interp1(depth,b_s(i,:,t),b_h(i).*b_sc(i,:));
        
        depth_t = depth;
        %depth_t(1) = depth(1) - s_el(i,t);
        b_u2_tem = ((b_u(i,1:end-1,t)+b_u(i,2:end,t))./2).*[(depth_t(2:end)-depth_t(1:end-1))]';
        b_u2(i,t) = sum(b_u2_tem(~isnan(b_u2_tem)));
     
        b_v2_tem = ((b_v(i,1:end-1,t)+b_v(i,2:end,t))./2).*[(depth_t(2:end)-depth_t(1:end-1))]';
        b_v2(i,t) = sum(b_v2_tem(~isnan(b_v2_tem)));       
    end
    
%-------------------------interpolate to s-level---------------------------
end
%-------------------------integrate to 2d----------------------------------

%------------------interpolate to out date---------------------------------

for i=1:b_num
    for n_l = 1:n_layer
        b_temp3(i,n_l,:) = interp1(data_day,squeeze(b_temp2(i,n_l,:)),date_out3d,'linear','extrap');
        b_s3(i,n_l,:) = interp1(data_day,squeeze(b_s2(i,n_l,:)),date_out3d,'linear','extrap');
    end
    b_u3(i,:) = interp1(data_day,squeeze(b_u2(i,:)),date_out2d,'linear','extrap');
    b_v3(i,:) = interp1(data_day,squeeze(b_v2(i,:)),date_out2d,'linear','extrap');
    b_el3(i,:) = interp1(data_day,b_el(i,:),date_out2d,'linear','extrap');
end
%------------------interpolate to out date---------------------------------       

%--------------------------rotate u,v to x,e-------------------------------


s_el3(:,1,:) = b_el3(sb_range,:);
w_el3(:,1,:) = b_el3(wb_range,:);
e_el3(:,1,:) = b_el3(eb_range,:);

s_el3(:,2,:) = b_el3(sb_range2,:);
w_el3(:,2,:) = b_el3(wb_range2,:);
e_el3(:,2,:) = b_el3(eb_range2,:);

s_el3(isnan(s_el3)) = 0;
w_el3(isnan(w_el3)) = 0;
e_el3(isnan(e_el3)) = 0;

s_s3(:,:,:) = b_s3(sb_range,:,:);
w_s3(:,:,:) = b_s3(wb_range,:,:);
e_s3(:,:,:) = b_s3(eb_range2,:,:);

s_temp3(:,:,:) = b_temp3(sb_range,:,:);
w_temp3(:,:,:) = b_temp3(wb_range,:,:);
e_temp3(:,:,:) = b_temp3(eb_range2,:,:);

s_u3(:,1,:) = b_u3(sb_range,:);
w_u3(:,1,:) = b_u3(wb_range,:);
e_u3(:,1,:) = b_u3(eb_range,:);

s_u3(:,2,:) = b_u3(sb_range2,:);
w_u3(:,2,:) = b_u3(wb_range2,:);
e_u3(:,2,:) = b_u3(eb_range2,:);

s_v3(:,1,:) = b_v3(sb_range,:);
w_v3(:,1,:) = b_v3(wb_range,:);
e_v3(:,1,:) = b_v3(eb_range,:);

s_v3(:,2,:) = b_v3(sb_range2,:);
w_v3(:,2,:) = b_v3(wb_range2,:);
e_v3(:,2,:) = b_v3(eb_range2,:);

%rotate u,v to x,e

[a,c,b] = size(s_u3);

%This suggests we need to process tidal boundary first
s_el_tide_s = load('../../tide_raw/s_el.mat');
w_el_tide_s = load('../../tide_raw/w_el.mat');
e_el_tide_s = load('../../tide_raw/e_el.mat');
s_el_tide = s_el_tide_s.s_el;
w_el_tide = w_el_tide_s.w_el;
e_el_tide = e_el_tide_s.e_el;



s_el_a = s_el_tide + squeeze(s_el3(:,1,:));
w_el_a = w_el_tide + squeeze(w_el3(:,1,:));
e_el_a = e_el_tide + squeeze(e_el3(:,2,:));

s_lon_u = rho2u_2d_mw(s_lon(:,1));
s_lon_v = rho2v_2d_mw(s_lon);
w_lon_u = rho2u_2d_mw(w_lon);
w_lon_v = rho2v_2d_mw(w_lon(1,:));
e_lon_u = rho2u_2d_mw(e_lon);
e_lon_v = rho2v_2d_mw(e_lon(2,:));

s_lat_u = rho2u_2d_mw(s_lat(:,1));
s_lat_v = rho2v_2d_mw(s_lat);
w_lat_u = rho2u_2d_mw(w_lat);
w_lat_v = rho2v_2d_mw(w_lat(1,:));
e_lat_u = rho2u_2d_mw(e_lat);
e_lat_v = rho2v_2d_mw(e_lat(2,:));

h_u = rho2u_2d_mw(h);
h_v = rho2v_2d_mw(h);

for i= 1:b
    s_el_a_u(:,i) = barnes(s_lon(:,1),s_lat(:,1),s_el_a(:,i),s_lon_u,s_lat_u,1,1,3);
    w_el_a_u(:,i) = barnes(w_lon(1,:),w_lat(1,:),w_el_a(:,i),w_lon_u,w_lat_u,1,1,3);
    e_el_a_u(:,i) = barnes(e_lon(2,:),e_lat(2,:),e_el_a(:,i),e_lon_u,e_lat_u,1,1,3);

    s_el_a_v(:,i) = barnes(s_lon(:,1),s_lat(:,1),s_el_a(:,i),s_lon_v,s_lat_v,1,1,3);
    w_el_a_v(:,i) = barnes(w_lon(1,:),w_lat(1,:),w_el_a(:,i),w_lon_v,w_lat_v,1,1,3);
    e_el_a_v(:,i) = barnes(e_lon(2,:),e_lat(2,:),e_el_a(:,i),e_lon_v,e_lat_v,1,1,3);
    
    s_u4(:,i) = rho2u_2d_mw(s_u3(:,1,i));
    w_u4(:,i) = rho2u_2d_mw(w_u3(:,:,i)');
    e_u4(:,i) = rho2u_2d_mw(e_u3(:,:,i)');
    
    s_v4(:,i) = rho2v_2d_mw(s_v3(:,:,i));
    w_v4(:,i) = rho2v_2d_mw(w_v3(:,1,i)');
    e_v4(:,i) = rho2v_2d_mw(e_v3(:,2,i)');
    
    s_uv3(:,i) = rho2u_2d_mw(s_v3(:,1,i));
    w_uv3(:,i) = rho2u_2d_mw(w_v3(:,:,i)');
    e_uv3(:,i) = rho2u_2d_mw(e_v3(:,:,i)');
    
    s_vu3(:,i) = rho2v_2d_mw(s_u3(:,:,i));
    w_vu3(:,i) = rho2v_2d_mw(w_u3(:,1,i)');
    e_vu3(:,i) = rho2v_2d_mw(e_u3(:,2,i)');  
    
    s_u4(:,i) = s_u4(:,i)./(h_u(:,1)+s_el_a_u(:,i));
    s_v4(:,i) = s_v4(:,i)./(h_v(:,1)+s_el_a_v(:,i));
    s_uv3(:,i) = s_uv3(:,i)./(h_u(:,1)+s_el_a_u(:,i));
    s_vu3(:,i) = s_vu3(:,i)./(h_v(:,1)+s_el_a_v(:,i));
    
    e_u4(:,i) = e_u4(:,i)./(h_u(end,:)'+e_el_a_u(:,i));
    e_v4(:,i) = e_v4(:,i)./(h_v(end,:)'+e_el_a_v(:,i));
    e_uv3(:,i) = e_uv3(:,i)./(h_u(end,:)'+e_el_a_u(:,i));
    e_vu3(:,i) = e_vu3(:,i)./(h_v(end,:)'+e_el_a_v(:,i));

    w_u4(:,i) = w_u4(:,i)./(h_u(1,:)'+w_el_a_u(:,i));
    w_v4(:,i) = w_v4(:,i)./(h_v(1,:)'+w_el_a_v(:,i));
    w_uv3(:,i) = w_uv3(:,i)./(h_u(1,:)'+w_el_a_u(:,i));
    w_vu3(:,i) = w_vu3(:,i)./(h_v(1,:)'+w_el_a_v(:,i));
end

s_angle_u = rho2u_2d_mw(s_angle(:,1));
w_angle_u = rho2u_2d_mw(w_angle);
e_angle_u = rho2u_2d_mw(e_angle);

s_angle_v = rho2v_2d_mw(s_angle);
w_angle_v = rho2v_2d_mw(w_angle(1,:));
e_angle_v = rho2v_2d_mw(e_angle(2,:));

s_mask_u = rho2u_2d_mw(s_mask(:,1));
w_mask_u = rho2u_2d_mw(w_mask);
e_mask_u = rho2u_2d_mw(e_mask);

s_mask_v = rho2v_2d_mw(s_mask);
w_mask_v = rho2v_2d_mw(w_mask(1,:));
e_mask_v = rho2v_2d_mw(e_mask(2,:));

s_el4(:,:) = squeeze(s_el3(:,1,:));
w_el4(:,:) = squeeze(w_el3(:,1,:));
e_el4(:,:) = squeeze(e_el3(:,1,:));

for t=1:b
    s_u5(:,t) = cos(s_angle_u).*s_u4(:,t)+sin(s_angle_u).*s_uv3(:,t);
    s_v5(:,t) = -sin(s_angle_v).*s_vu3(:,t)+cos(s_angle_v).*s_v4(:,t);
    e_u5(:,t) = cos(e_angle_u').*e_u4(:,t)+sin(e_angle_u').*e_uv3(:,t);
    e_v5(:,t) = -sin(e_angle_v').*e_vu3(:,t)+cos(e_angle_v').*e_v4(:,t);
    w_u5(:,t) = cos(w_angle_u').*w_u4(:,t)+sin(w_angle_u').*w_uv3(:,t);
    w_v5(:,t) = -sin(w_angle_v').*w_vu3(:,t)+cos(w_angle_v').*w_v4(:,t);
          
    s_u5(s_mask_u==0,t) = 0;
    s_v5(s_mask_v==0,t) = 0;

    e_u5(e_mask_u==0,t) = 0;
    e_v5(e_mask_v==0,t) = 0;

    w_u5(w_mask_u==0,t) = 0;
    w_v5(w_mask_v==0,t) = 0;
end
%--------------------------rotate u,v to x,e-------------------------------

%Output

s_el_out = s_el4;
s_s_out = s_s3;
s_temp_out = s_temp3;
s_u_out = s_u5;
s_v_out = s_v5;

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

save('s_el_out.mat','s_el_out');           
save('s_s_out.mat','s_s_out');
save('s_temp_out.mat','s_temp_out');
save('s_u_out.mat','s_u_out');
save('s_v_out.mat','s_v_out');

save('w_el_out.mat','w_el_out');           
save('w_s_out.mat','w_s_out');
save('w_temp_out.mat','w_temp_out');
save('w_u_out.mat','w_u_out');
save('w_v_out.mat','w_v_out');

save('e_el_out.mat','e_el_out');           
save('e_s_out.mat','e_s_out');
save('e_temp_out.mat','e_temp_out');
save('e_u_out.mat','e_u_out');
save('e_v_out.mat','e_v_out');


%Figures for checking
figure;
quiver(s_lon_u,s_lat_u,s_u4(:,1),s_v4(1:end-1,1),0,'r');
hold on;
quiver(s_lon_u,s_lat_u,s_u_out(:,1,1),s_v_out(1:end-1,1,1),0);
title('Soutern Boundary Velocity');
legend('Before rotation','After rotation');

figure;
quiver(w_lon_v',w_lat_v',w_u4(1:end-1,1),w_v4(:,1),0,'k');
hold on;
quiver(w_lon_v',w_lat_v',w_u_out(1:end-1,1,1),w_v_out(:,1,1),0,'g');
title('Western Boundary Velocity');
legend('Before rotation','After rotation');

figure;
quiver(e_lon_v',e_lat_v',e_u4(1:end-1,1),e_v4(:,1),0,'r');
hold on;
quiver(e_lon_v',e_lat_v',e_u_out(1:end-1,1,1),e_v_out(:,1,1),0);
title('Eastern Boundary Velocity');
legend('Before rotation','After rotation');

figure;
subplot(1,3,1);
contourf(w_temp3(:,:,1));
title('Western Boundary Temp.');
subplot(1,3,2);
contourf(e_temp3(:,:,1));
title('Eastern Boundary Temp.');
subplot(1,3,3);
contourf(s_temp3(:,:,1));
title('Soutern Boundary Temp.');

figure;
subplot(1,3,1);
contourf(w_s3(:,:,20));
title('Western Boundary Sal.');
subplot(1,3,2);
contourf(e_s3(:,:,20));
title('Eastern Boundary Sal.');
subplot(1,3,3);
contourf(s_s3(:,:,20));
title('Soutern Boundary Sal.');




