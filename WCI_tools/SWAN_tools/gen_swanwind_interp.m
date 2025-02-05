clear all;
grd = 'ROMS_WFS_new.nc';
lon = ncread(grd,'lon_rho');
lat = ncread(grd,'lat_rho');
lat = lat';
lon = lon';
year = 2022;
fid3 = fopen(['WFS_swan_',num2str(year),'_wind_no_IAN.dat'],'wt');
outname1 = '../No_Hurricane_File/WFS_2022_Uwind1_no_IAN.nc';
outname2 = '../No_Hurricane_File/WFS_2022_Vwind1_no_IAN.nc';

u_out = ncread(outname1,'Uwind');
v_out = ncread(outname2,'Vwind');

[~,~,n] = size(u_out);

u_final(:,:,1:n) = u_out;
v_final(:,:,1:n) = v_out;

origin_date = double(ncread(outname1,'wind_time'))+datenum(year,1,1);
out_date = datenum(year,1,1,0,0,0):6/24:datenum(year,12,31,24,0,0);

u_final2 = tv3dinterpt(u_final,origin_date,out_date,'linear');
v_final2 = tv3dinterpt(v_final,origin_date,out_date,'linear');

[r,c] = size(lat);

for i=1:size(u_final2,3)
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