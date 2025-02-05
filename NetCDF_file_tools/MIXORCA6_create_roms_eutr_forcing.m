clear all;close all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\NC_file_generation');
grd = '../Model_grid/ROMS_WFS_new.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
mask = ncread(grd,'mask_rho');
gn = struct('lon_rho',lon);
gn.N = length(ncread(grd,'Cs_r'));
fn1 = 'WFS_2005_2006_Uwind1_bio_mixo.nc';
fn2 = 'WFS_2005_2006_Uwind2_bio_mixo.nc';
fn3 = 'WFS_2005_2006_Vwind1_bio_mixo.nc';
fn4 = 'WFS_2005_2006_Vwind2_bio_mixo.nc';
fn5 = 'WFS_2005_2006_swrad1_bio_mixo.nc';
fn6 = 'WFS_2005_2006_swrad2_bio_mixo.nc';
fn7 = 'WFS_2005_2006_PAR_bio_mixo.nc';

[r,c] = size(lon);

time_ref = datenum(2005,6,1,0,0,0);

wind_intens = 1.;

origin_date = datenum(2005,1,1,0,0,0):6/24:datenum(2006,12,31,24,0,0);

out_date1 = datenum(2005,6,1,0,0,0):6/24:datenum(2006,2,1,1,0,0);
out_date2 = [];

time1 = out_date1-time_ref;
time2 = out_date2-time_ref;

year = datevec(out_date1(1));
year = year(1);

outname1 = strcat('atm_forcing_',num2str(year),'.mat');
outname2 = strcat('atm_forcing_',num2str(year+1),'.mat');
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


u1 = load(strcat('../Atmosphere_forcing_preprocessing/',outname1),'u_out');
u2 = load(strcat('../Atmosphere_forcing_preprocessing/',outname2),'u_out');
u_out = u1.u_out;
u_out(:,:,end:end+size(u2.u_out,3)-1) = u2.u_out;

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
clear u1 u2 u_out u_out u_out2 Uwind 



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

v1 = load(strcat('../Atmosphere_forcing_preprocessing/',outname1),'v_out');
v2 = load(strcat('../Atmosphere_forcing_preprocessing/',outname2),'v_out');
v_out = v1.v_out;
v_out(:,:,end:end+size(v2.v_out,3)-1) = v2.v_out;

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
clear v1 v2 v_out v_out v_out2 Vwind 


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


ssr_d1 = load(strcat('../Atmosphere_forcing_preprocessing/',outname1),'ssr_d_out');
ssr_d2 = load(strcat('../Atmosphere_forcing_preprocessing/',outname2),'ssr_d_out');
ssr_d_out = ssr_d1.ssr_d_out;
ssr_d_out(:,:,end:end+size(ssr_d2.ssr_d_out,3)-1) = ssr_d2.ssr_d_out;

ssr_u1 = load(strcat('../Atmosphere_forcing_preprocessing/',outname1),'ssr_u_out');
ssr_u2 = load(strcat('../Atmosphere_forcing_preprocessing/',outname2),'ssr_u_out');
ssr_u_out = ssr_u1.ssr_u_out;
ssr_u_out(:,:,end:end+size(ssr_u2.ssr_u_out,3)-1) = ssr_u2.ssr_u_out;

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
clear ssr_d1 ssr_d2 ssr_u1 ssr_u2 swrad ssr_out 

%PAR
par1 = load(['../PAR/','PAR_',num2str(year),'.mat']);
par2 = load(['../PAR/','PAR_',num2str(year+1),'.mat']);
PAR = par1.PAR;
PAR(:,:,end:end+size(par2.PAR,3)-1) = par2.PAR;
tday1 = par1.tday+datenum(year,1,1);
tday2 = par2.tday+datenum(year+1,1,1);
tday = [tday1(:,1:end-1),tday2];

pos = find(tday-time_ref>=0);
tday2 = tday(pos)-time_ref;
PAR = PAR(:,:,pos);

create_roms_forcings(lon,lat,tday2,fn7,'PAR');