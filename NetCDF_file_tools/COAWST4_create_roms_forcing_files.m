clear all;close all;
bay_mask = ncread('../Model_grid/ROMS_WFS_10river_grid_bay_mask.nc','mask_rho');
grd = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd,'Cs_r'));
dep = ncread(grd,'h');
year = 2021;
fn1 = ['WFS_',num2str(year),'_Uwind1.nc'];
fn2 = ['WFS_',num2str(year),'_Uwind2.nc'];
fn3 = ['WFS_',num2str(year),'_Vwind1.nc'];
fn4 = ['WFS_',num2str(year),'_Vwind2.nc'];
fn5 = ['WFS_',num2str(year),'_Pair1.nc'];
fn6 = ['WFS_',num2str(year),'_Pair2.nc'];
fn7 = ['WFS_',num2str(year),'_Tair1.nc'];
fn8 = ['WFS_',num2str(year),'_Tair2.nc'];
fn9 = ['WFS_',num2str(year),'_Qair1.nc'];
fn10 = ['WFS_',num2str(year),'_Qair2.nc'];
fn11 = ['WFS_',num2str(year),'_swrad1.nc'];
fn12 = ['WFS_',num2str(year),'_swrad2.nc'];
fn13 = ['WFS_',num2str(year),'_lwrad1.nc'];
fn14 = ['WFS_',num2str(year),'_lwrad2.nc'];
fn15 = ['WFS_',num2str(year),'_rain1.nc'];
fn16 = ['WFS_',num2str(year),'_rain2.nc'];

scale_depth = get_scale_depth(dep,bay_mask,mask,0.82,0.9956);
%wind_intens = 1.3*scale_depth;
wind_intens = 1.;

origin_date = datenum(year,1,1,0,0,0):6/24:datenum(year,12,31,24,0,0);

out_date1 = datenum(year,1,1,0,0,0):6/24:datenum(year,12,31,24,0,0);
out_date2 = [];

time1 = out_date1-out_date1(1);
time2 = out_date2-out_date1(1);

outname = strcat('atm_forcing_',num2str(year),'.mat');

%UWIND
clear out_i1 out_i2
for i=1:length(out_date1)
    delta = abs(origin_date(:)-out_date1(i));
    out_i1(i) = find(delta==min(delta));
end

for i=1:length(out_date2)
    delta = abs(origin_date(:)-out_date2(i));
    out_i2(i) = find(delta==min(delta));
end


load(strcat('../Atmosphere_forcing_preprocessing/',outname),'u_out');
%u_out2 = u_out;
for ii=1:size(u_out,3)
    u_out2(:,:,ii) = u_out(:,:,ii).*wind_intens;
end

Uwind = u_out2(:,:,out_i1);
if(sum(sum(sum(isnan(Uwind))))>0)
    pause;
end
create_roms_forcings(lon,lat,time1,fn1,'Uwind');

if(~isempty(time2))
Uwind = u_out2(:,:,out_i2);
if(sum(sum(sum(isnan(Uwind))))>0)
    pause;
end
create_roms_forcings(lon,lat,time2,fn2,'Uwind');
end
clear u_out u_out u_out2 Uwind 



%VWIND
clear out_i1 out_i2
for i=1:length(out_date1)
    delta = abs(origin_date(:)-out_date1(i));
    out_i1(i) = find(delta==min(delta));
end

for i=1:length(out_date2)
    delta = abs(origin_date(:)-out_date2(i));
    out_i2(i) = find(delta==min(delta));
end

load(strcat('../Atmosphere_forcing_preprocessing/',outname),'v_out');
%v_out2 = v_out;

for ii=1:size(v_out,3)
    v_out2(:,:,ii) = v_out(:,:,ii).*wind_intens;
end

Vwind = v_out2(:,:,out_i1);
if(sum(sum(sum(isnan(Vwind))))>0)
    pause;
end
create_roms_forcings(lon,lat,time1,fn3,'Vwind');

if(~isempty(time2))
Vwind = v_out2(:,:,out_i2);
if(sum(sum(sum(isnan(Vwind))))>0)
    pause;
end
create_roms_forcings(lon,lat,time2,fn4,'Vwind');
end
clear v_out v_out v_out2 Vwind 


%PAIR
clear out_i1 out_i2
for i=1:length(out_date1)
    delta = abs(origin_date(:)-out_date1(i));
    out_i1(i) = find(delta==min(delta));
end

for i=1:length(out_date2)
    delta = abs(origin_date(:)-out_date2(i));
    out_i2(i) = find(delta==min(delta));
end

load(strcat('../Atmosphere_forcing_preprocessing/',outname),'slp_out');
slp_out = slp_out/100;
slp_out2 = slp_out;

Pair = slp_out2(:,:,out_i1);
if(sum(sum(sum(isnan(Pair))))>0)
    pause;
end
figure; plot(squeeze(Pair(100,100,:)));
create_roms_forcings(lon,lat,time1,fn5,'Pair');

if(~isempty(time2))
Pair = slp_out2(:,:,out_i2);
if(sum(sum(sum(isnan(Pair))))>0)
    pause;
end
create_roms_forcings(lon,lat,time2,fn6,'Pair');
figure; plot(squeeze(Pair(100,100,:)));
end
clear slp_out slp_out slp_out2 Pair

%TAIR
clear out_i1 out_i2
for i=1:length(out_date1)
    delta = abs(origin_date(:)-out_date1(i));
    out_i1(i) = find(delta==min(delta));
end

for i=1:length(out_date2)
    delta = abs(origin_date(:)-out_date2(i));
    out_i2(i) = find(delta==min(delta));
end

load(strcat('../Atmosphere_forcing_preprocessing/',outname),'t_out');
t_out = t_out-273.15;
Tair = t_out(:,:,out_i1);
if(sum(sum(sum(isnan(Tair))))>0)
    pause;
end
figure; plot(squeeze(Tair(100,100,:)));
create_roms_forcings(lon,lat,time1,fn7,'Tair');

if(~isempty(time2))
Tair = t_out(:,:,out_i2);
if(sum(sum(sum(isnan(Tair))))>0)
    pause;
end
figure; plot(squeeze(Tair(100,100,:)));
create_roms_forcings(lon,lat,time2,fn8,'Tair');
end
clear t_out Tair t_out

%HUMIDITY
clear out_i1 out_i2
for i=1:length(out_date1)
    delta = abs(origin_date(:)-out_date1(i));
    out_i1(i) = find(delta==min(delta));
end

for i=1:length(out_date2)
    delta = abs(origin_date(:)-out_date2(i));
    out_i2(i) = find(delta==min(delta));
end

load(strcat('../Atmosphere_forcing_preprocessing/',outname),'hum_out');
Qair = hum_out(:,:,out_i1);
if(sum(sum(sum(isnan(Qair))))>0)
    pause;
end
figure; plot(squeeze(Qair(100,100,:)));
create_roms_forcings(lon,lat,time1,fn9,'Qair');

if(~isempty(time2))
Qair = hum_out(:,:,out_i2);
if(sum(sum(sum(isnan(Qair))))>0)
    pause;
end
figure; plot(squeeze(Qair(100,100,:)));
create_roms_forcings(lon,lat,time2,fn10,'Qair');
end
clear Qair hum_out 

%SSR
clear out_i1 out_i2
for i=1:length(out_date1)
    delta = abs(origin_date(:)-out_date1(i));
    out_i1(i) = find(delta==min(delta));
end

for i=1:length(out_date2)
    delta = abs(origin_date(:)-out_date2(i));
    out_i2(i) = find(delta==min(delta));
end

load(strcat('../Atmosphere_forcing_preprocessing/',outname),'ssr_d_out');
load(strcat('../Atmosphere_forcing_preprocessing/',outname),'ssr_u_out');
ssr_out = ssr_d_out - ssr_u_out;
%%%%%%%%%%%ssr correction%%%%%%%%%%%
%ssr_out = ssr_out .* 0.85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
swrad = ssr_out(:,:,out_i1);
if(sum(sum(sum(isnan(swrad))))>0)
    pause;
end
figure; plot(squeeze(swrad(100,100,:)));
create_roms_forcings(lon,lat,time1,fn11,'swrad');

if(~isempty(time2))
swrad = ssr_out(:,:,out_i2);
if(sum(sum(sum(isnan(swrad))))>0)
    pause;
end
figure; plot(squeeze(swrad(100,100,:)));
create_roms_forcings(lon,lat,time2,fn12,'swrad');
end
clear swrad ssr_out 

%SLR_D
clear out_i1 out_i2
for i=1:length(out_date1)
    delta = abs(origin_date(:)-out_date1(i));
    out_i1(i) = find(delta==min(delta));
end

for i=1:length(out_date2)
    delta = abs(origin_date(:)-out_date2(i));
    out_i2(i) = find(delta==min(delta));
end

load(strcat('../Atmosphere_forcing_preprocessing/',outname),'slr_d_out');
%%%%%%%%%%%slr correction%%%%%%%%%%%
%slr_d_out = slr_d_out .* 0.85;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lwrad = slr_d_out(:,:,out_i1);
lwrad_down = lwrad;
if(sum(sum(sum(isnan(lwrad_down))))>0)
    pause;
end
figure; plot(squeeze(lwrad_down(100,100,:)));
create_roms_forcings(lon,lat,time1,fn13,'lwrad_down');

if(~isempty(time2))
lwrad = slr_d_out(:,:,out_i2);
lwrad_down = lwrad;
if(sum(sum(sum(isnan(lwrad_down))))>0)
    pause;
end
figure; plot(squeeze(lwrad_down(100,100,:)));
create_roms_forcings(lon,lat,time2,fn14,'lwrad_down');
end
clear lwrad_down lwrad slr_d_out

clear all;
% RAIN
% clear out_i1 out_i2
% for i=1:length(out_date1)
%     delta = abs(origin_date(:)-out_date1(i));
%     out_i1(i) = find(delta==min(delta));
% end
% 
% for i=1:length(out_date2)
%     delta = abs(origin_date(:)-out_date2(i));
%     out_i2(i) = find(delta==min(delta));
% end
% 
% load('../Atmosphere_forcing_preprocessing/rain_out.mat');
% rain = rain_out(:,:,out_i1);
% if(sum(sum(sum(isnan(rain))))>0)
%     pause;
% end
% figure; plot(squeeze(rain(100,100,:)));
% create_roms_forcings(lon,lat,time1,fn15,'rain');
% 
% if(~isempty(time2))
% rain = rain_out(:,:,out_i2);
% if(sum(sum(sum(isnan(rain))))>0)
%     pause;
% end
% figure; plot(squeeze(rain(100,100,:)));
% create_roms_forcings(lon,lat,time2,fn16,'rain');
% end
% clear rain rain_out
