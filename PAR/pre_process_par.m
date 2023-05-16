clear all;close all;

load('C:\Users\cheny\Desktop\EcoHAB\PAR\erdMH1par01day_gom_2003_2013.mat')
time = erdMH1par01day.time/3600/24+datenum(1970,1,1);
par_lon = double(erdMH1par01day.longitude);
par_lat = double(erdMH1par01day.latitude);
[plon,plat] = meshgrid(par_lon,par_lat);
par0 = double(erdMH1par01day.par);

load('C:\Users\cheny\Desktop\EcoHAB\PAR\erdMH1par01day_gom_2013_2021.mat')
time = [time; erdMH1par01day.time/3600/24+datenum(1970,1,1)];
par0(end+1:length(time),:,:) = double(erdMH1par01day.par);

par0 = par0*1e6/24/3600 ;%Einsteins m-2 d-1 -> µmoles m-2 s-1
tvec = datevec(time);

[r,c] = size(plon);

par1 = zeros(length(time),r,c);
for i=1:r
    for j=1:c
        y = par0(:,i,j);
        if(sum(~isnan(y))>0)
            par1(:,i,j) = interp1(time(~isnan(y)),y(~isnan(y)),time,'linear');
            y = par1(:,i,j);
            par1(:,i,j) = interp1(time(~isnan(y)),y(~isnan(y)),time,'nearest','extrap');
        else
            par1(:,i,j) = NaN;
        end
    end
end

for year=2003:2021
    pos = find(tvec(:,1)==year);
    par = par1(pos,:,:);
    save(['./PAR_raw/PAR_raw_',num2str(year),'.mat'],"time","par","plat","plon","-v7.3");
end
