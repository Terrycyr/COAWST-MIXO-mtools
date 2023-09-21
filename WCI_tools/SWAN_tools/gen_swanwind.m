clear all;
lat = ncread('./pre_2017_grd.nc','lat_rho');
lon = ncread('./pre_2017_grd.nc','lon_rho');
lat = lat';
lon = lon';
year = 2017;
days = 62 ;
fid3 = fopen(['PRE_',num2str(year),'_wind.dat'],'wt');

load('./wind_raw/u_final.mat');
load('./wind_raw/v_final.mat');
load('./wind_raw/u_out.mat');
load('./wind_raw/v_out.mat');

[r,c,n] = size(u_out);
[r1,c1,n1] = size(u_final);

u_final2(:,:,1:n) = u_out;
u_final2(:,:,n+1:n1+n) = u_final;

v_final2(:,:,1:n) = v_out;
v_final2(:,:,n+1:n1+n) = v_final;

%original 7.1 to 9.30
origin_date = datenum(2017,7,1,0,0,0):6/24:datenum(2017,9,30,24,0,0);
out_date = datenum(2017,7,20,0,0,0):6/24:datenum(2017,9,19,24,0,0);

for i=1:length(out_date)
    out_i(i) = find(origin_date(:)==out_date(i));
end

u_final2 = u_final2(:,:,out_i);
v_final2 = v_final2(:,:,out_i);

[r,c] = size(lat);

for i=1:days*4+1
    i
    for ii=1:r
        fprintf(fid3,'%10.3f',u_final2(:,ii,i));
        fprintf(fid3,'\n');
    end

    for ii=1:r
        fprintf(fid3,'%10.3f',v_final2(:,ii,i));
        fprintf(fid3,'\n');
    end
   % fwrite(fid4,eps,'single');
   % fwrite(fid4,u_final,'single');
   % fwrite(fid4,v_final,'single');
   % fwrite(fid4,eps,'single');
end

fclose(fid3);