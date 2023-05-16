clear all;
fn = './trans_list.dat';
fid = fopen(fn,'wt+');
year = 2001;

for i=1:365
    if(i==1)
        fprintf(fid,'     CLMNAME == ROMS_TRANS/WFS_%s_his_%s.nc |\n',num2str(year),sprintf('%04d',i));
    elseif(i==365)
        fprintf(fid,'                ROMS_TRANS/WFS_%s_his_%s.nc ',num2str(year),sprintf('%04d',i));
    else
        fprintf(fid,'                ROMS_TRANS/WFS_%s_his_%s.nc |\n',num2str(year),sprintf('%04d',i));
    end
end

fclose(fid);