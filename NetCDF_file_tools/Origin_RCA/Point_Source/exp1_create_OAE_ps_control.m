clear all; close all;
addpath(path,'C:\Users\cheny\Desktop\matlab_tools\fill_coastlines');
grd_name = 'cpb_grd_80x120_wf32.nc';

filename = "C:\Users\cheny\Desktop\matlab_tools\m_map\private\gshhs_h.b";
S = gshhs(filename{1},[26.2 27.3],[-82.9 -81.9]);
clat = [S.Lat];
clon = [S.Lon];

lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');
gn = struct('lon_rho',lon);

ps_fname = 'CPB_WWTP_ps_oae_control.inp';

SYSNAME = {"SAL","PHYT1","PHYT2","PHYT3","RPOP","LPOP","RDOP","LDOP","PO4T","RPON"...
    ,"LPON","RDON","LDON","NH4T","NO23","BSI","SIT","RPOC","LPOC","RDOC"...
    ,"LDOC","EXDOC","REPOC","REDOC","O2EQ","DO"}; %#ok<*CLARRSTR>

unit = {'kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day'...
    ,'kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','kg/day'...
    ,'kg/day','kg/day','kg/day','kg/day','kg/day','kg/day','mol/day','mol/day','kg/day','kg/day'};

year = 2010;

IPSOPT = 2;
IPSPWLOPT = 1;
N = 20;
lb2kg = 0.453592;
MG2M3 = 3785.41178;
MG2L = 3785411.78;
wdensity = 1000; %kg/m3
wdensity2 = 1; %kg/L
k1k2_m = 9;

%------------------------------Constant Data-------------------------------
temp = 20;
TA_const = 2120; %umol/kg
DIC_const = 3628; %umol/kg
NH3_const = 0.47; %N mg/L
NO23_const = 1.19; %N mg/L
PO4_const = 0.; %P mg/L
SI_const = 0.; %SI mg/L
sal_const = 0.1;
pH = cal_pH(TA_const,DIC_const,sal_const,temp,temp,0,0,SI_const/28/wdensity2*1000,PO4_const/31/wdensity2*1000,k1k2_m);
%--------------------------------------------------------------------------

load('WWTP_flow.mat');
station_id = 1:size(location,1);

figure;
%contourf(lon,lat,mask);
%hold on;
fill_coastline(clon,clat);
axis image;
axis([-77.46 -75.05 36.23 39.56]);
hold on;
for i=1:length(station_ij)
    scatter(lon(station_ij(i,1),station_ij(i,2)),lat(station_ij(i,1),station_ij(i,2)),40,'green','filled');
    hold on;
    text(lon(station_ij(i,1),station_ij(i,2)),lat(station_ij(i,1),station_ij(i,2)),location_name{i},'color','r','fontweight','bold');
end
set(gcf,'color','w');
xlabel('Longitude');
ylabel('Latitude');

% Fraction in each sigma layer
ZFRACPS = [0.5 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%UNIT OF TBREAK/ps_t: SECS, MINS, HRS, DAYS
TWARPPS = 'DAYS';

ps_t = [datenum(year-1,12,1),datenum(year,1:12,1),datenum(year+1,1,1)] - datenum(1983,1,1);
ps = zeros(size(location,1),length(ps_t),length(SYSNAME));

for i=1:length(ps_t)
    tmp = datevec(ps_t(i)+datenum(1983,1,1));
    mon = tmp(2);
    ps(:,i,14) = NH3_const.*FLOW_mean*MG2L/1e6;
    ps(:,i,15) = NO23_const.*FLOW_mean*MG2L/1e6;
    ps(:,i,9) = PO4_const.*FLOW_mean*MG2L/1e6;
    ps(:,i,23) = TA_const.*FLOW_mean*MG2M3*wdensity/1e6;  %TA umol kg-1 => mol day-1
    ps(:,i,24) = DIC_const.*FLOW_mean*MG2M3*wdensity/1e6; %DIC umol kg-1 => mol day-1
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
    fprintf(fid,'%5i%5i%5i%50s\n',IPS(i,1),IPS(i,2),IPS(i,3),location_name{i});
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