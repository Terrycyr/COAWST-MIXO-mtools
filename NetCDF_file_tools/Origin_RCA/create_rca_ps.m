clear all; close all;

grd_name = 'cpb_grd_80x120_wf32.nc';

lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');
gn = struct('lon_rho',lon);

ps_fname = 'CPB_WWTP_ps_test.inp';

SYSNAME = {"SAL","PHYT1","PHYT2","PHYT3","RPOP","LPOP","RDOP","LDOP","PO4T","RPON"...
    ,"LPON","RDON","LDON","NH4T","NO23","BSI","SIT","RPOC","LPOC","RDOC"...
    ,"LDOC","EXDOC","REPOC","REDOC","O2EQ","DO"}; %#ok<*CLARRSTR>

unit = {'kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day'...
    ,'kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day'...
    ,'kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','mol/day','mol/day','kg/day','kg/day'};

IPSOPT = 2;
IPSPWLOPT = 1;
N = 20;
lb2kg = 0.453592;
MG2M3 = 3785.41178;
wdensity = 1000;

location = {'HRSD-NANSEMOND';'Blue Plains';'Cambridge';'Back River';'PATAPSCO'};
station_id = 1:5;
station_ij = [54    7;...
               4   83;...
              78   79;...
              34  105;...
              30  105];    

figure;
contourf(lon,lat,mask);
hold on;
for i=1:length(station_ij)
    scatter(lon(station_ij(i,1),station_ij(i,2)),lat(station_ij(i,1),station_ij(i,2)),40,'green','filled');
    hold on;
    text(lon(station_ij(i,1),station_ij(i,2)),lat(station_ij(i,1),station_ij(i,2)),location{i});
end

% Fraction in each sigma layer
ZFRACPS = [0.5 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%UNIT OF TBREAK/ps_t: SECS, MINS, HRS, DAYS
TWARPPS = 'DAYS';

%unit of ps: kg/day !!!
[data_time,BOD5,DO,FLOW,NH3,NO23,PO4,TKN,TN,TON,TOP,TP,TSS] = get_ps('2000 Point Sources MD.DE.DC.VA.xlsx',location);

for i=1:length(location)
    FLOW_mon(i,:) = interp1(data_time{i}{3},FLOW{i},datenum(2000,1:12,1),'nearest','extrap');
    NH3_mon(i,:) = interp1(data_time{i}{4},NH3{i},datenum(2000,1:12,1),'nearest','extrap');
    NO23_mon(i,:) = interp1(data_time{i}{5},NO23{i},datenum(2000,1:12,1),'nearest','extrap');
    PO4_mon(i,:) = interp1(data_time{i}{6},PO4{i},datenum(2000,1:12,1),'nearest','extrap');
end


ps_t = [datenum(1999,12,1),datenum(2000,1:12,1),datenum(2001,1,1)] - datenum(1983,1,1);
ps = zeros(length(location),length(ps_t),length(SYSNAME));

for i=1:length(ps_t)
    tmp = datevec(ps_t(i)+datenum(1983,1,1));
    mon = tmp(2);
    ps(:,i,14) = NH3_mon(:,mon)*lb2kg;
    ps(:,i,15) = NO23_mon(:,mon)*lb2kg;
    ps(:,i,9) = PO4_mon(:,mon)*lb2kg;
    ps(:,i,23) = 3000.*FLOW_mon(:,mon)*MG2M3*wdensity/1e6;  %TA 3000 umol kg-1 => mol day-1
    ps(:,i,24) = 3000/0.93.*FLOW_mon(:,mon)*MG2M3*wdensity/1e6;  %DIC TA/0.93
end


%NOPS = all, SCALE = 1.0
IPS = [station_id',station_ij];
NOPS = size(IPS,1)*ones(1,length(SYSNAME));
SCALPS = ones(1,length(SYSNAME));



%--------------------------------------------------------------------------
%---------------------------- NO CHANGE BELOW -----------------------------
%--------------------------------------------------------------------------

fid =  fopen(ps_fname,'wt+');

fprintf(fid,'%s\n','         TV        PL (TV:1=const., 2=time var., PL:0=step-func, 1=piece. linear)');
fprintf(fid,'%10i%10i\n',IPSOPT, IPSPWLOPT);
fprintf(fid,'%s\n','    NO   IX   IY Location / Vertical Distribution Fractions');
for i=1:size(IPS,1)
    fprintf(fid,'%5i%5i%5i%50s\n',IPS(i,1),IPS(i,2),IPS(i,3),location{i});
    fprintf(fid,[repmat(' ',1,10),repmat('%6.3f',1,N),'\n'], ZFRACPS);
end
fprintf(fid,'%s\n','  -99           End of Loading Table Index');
fprintf(fid,'%s\n','              TBREAK   UNIT');


for t = 1:length(ps_t)
    for i = 1:size(SYSNAME,2)
        if(t==1)
            if(i==1)
                fprintf(fid,[repmat(' ',1,10),'%10i',repmat(' ',1,6),'%4s\n'], ps_t(t),TWARPPS);
            end
            fprintf(fid,'%s\n',append('No. of Discharge for system ',num2str(i),' Loads,  ',SYSNAME{i},' ',unit{i}));
            fprintf(fid,'%10i\n',NOPS(i));
            fprintf(fid,'%s\n','Scale Factor');
            fprintf(fid,'%10.3e\n',SCALPS(i));

            if(size(IPS,1)<7)
                fprintf(fid,[repmat(' ',1,10),repmat('%10i',1,size(IPS,1)),'\n'],IPS(:,1));
            else
                tmp = mod(size(IPS,1),7);
                for j = 1:(size(IPS,1)-tmp)/7
                    fprintf(fid,[repmat(' ',1,10),repmat('%10i',1,7),'\n'],IPS([1:7]+(j-1)*7,1));
                end
                fprintf(fid,[repmat(' ',1,10),repmat('%10i',1,tmp),'\n'],IPS(end-tmp+1:end,1));
            end  
        else
            if(i==1)
                fprintf(fid,[repmat(' ',1,10),'%10i',repmat(' ',1,6),'%6s\n'],ps_t(t),'TBREAK');
            end
        end

        if(size(IPS,1)<7)
            fprintf(fid,[repmat(' ',1,10),repmat('%12.3f',1,size(IPS,1)),'\n'],ps(:,t,i));
        else
            tmp = mod(size(IPS,1),7);
            for j = 1:(size(IPS,1)-tmp)/7
                fprintf(fid,[repmat(' ',1,10),repmat('%12.3f',1,7),'\n'],ps([1:7]+(j-1)*7,t,i));
            end
            fprintf(fid,[repmat(' ',1,10),repmat('%12.3f',1,tmp),'\n'],ps(end-tmp+1:end,t,i));
        end
    end
end
fclose(fid);