clear all;close all;
grd = '../Model_grid/ROMS_WFS_10river_grid_v11.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd,'Cs_r'));
fn1 = 'WFS_2001_Uwind1_bio.nc';
fn2 = 'WFS_2001_Uwind2_bio.nc';
fn3 = 'WFS_2001_Vwind1_bio.nc';
fn4 = 'WFS_2001_Vwind2_bio.nc';
fn5 = 'WFS_2001_swrad1_bio.nc';
fn6 = 'WFS_2001_swrad2_bio.nc';
fn7 = 'WFS_2001_PAR_bio.nc';

[r,c] = size(lon);

time_ref = datenum(2001,4,1,0,0,0);

wind_intens = 1.;

origin_date = datenum(2001,1,1,0,0,0):6/24:datenum(2001,12,31,24,0,0);

out_date1 = datenum(2001,1,1,0,0,0):6/24:datenum(2001,12,31,24,0,0);
out_date2 = [];

time1 = out_date1-time_ref;
time2 = out_date2-time_ref;

year = datevec(out_date1(1));
year = year(1);
outname = strcat('atm_forcing_',num2str(year),'.mat');

time_ref2 = time_ref - datenum(year,1,1);

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
create_roms_forcings(lon,lat,time1,fn5,'swrad');

if(~isempty(time2))
swrad = ssr_out(:,:,out_i2);
if(sum(sum(sum(isnan(swrad))))>0)
    pause;
end
figure; plot(squeeze(swrad(100,100,:)));
create_roms_forcings(lon,lat,time2,fn6,'swrad');
end
clear swrad ssr_out 

%PAR
load(['../PAR/','PAR_',num2str(year),'.mat']);
create_roms_forcings(lon,lat,tday-time_ref2,fn7,'PAR');
