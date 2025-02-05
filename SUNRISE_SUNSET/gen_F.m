clear all;

date_out = datenum(2002,1,1,0,0,0):1:datenum(2002,12,31,24,0,0);

for i=1:length(date_out)
    
    F(i) = sunrise(27.05,-82.75,0,0,datestr(date_out(i),'yyyy-mm-dd'));
    
end

save('F.mat','F');