clear all; close all;

grd_name = 'cpb_grd_80x120_wf32.nc';

lon = ncread(grd_name,'lon_rho');
lat = ncread(grd_name,'lat_rho');
mask = ncread(grd_name,'mask_rho');
h = ncread(grd_name,'h');
gn = struct('lon_rho',lon);

year = 2010;

IPSOPT = 2;
IPSPWLOPT = 1;
N = 20;
lb2kg = 0.453592;
MG2M3 = 3785.41178;
MG2L = 3785411.78;

location = {{'HRSD-VIP';...
            'HRSD-NANSEMOND';...
            'HRSD-ARMY BASE';...
            'HRSD-BOAT HARBOR';...
            'HRSD-JAMES RIVER';...
            'HRSD-WILLIAMSBURG';...
            'HRSD-CHESAPEAKE/ELIZABETH';...
            'HRSD-YORK'}
            {'Blue Plains';'WASHINGTON, D.C.  COMBINED SEWER OVERFLOW';...
            'NSWC-INDIAN HEAD';...
            'INDIAN HEAD';...
            'NSWC-INDIAN HEAD';...
            'MATTAWOMAN';...
            'LACKEY HIGH';...
            'DALE CITY #8';...
            'DALE CITY #1';...
            'H.L. MOONEY';...
            'ARLINGTON';...
            'ALEXANDRIA';...
            'NOMAN M. COLE JR. POLLUTION CONTROL PLANT';...
            'QUANTICO-MAINSIDE';...
            'AQUIA'}
            {'Back River';...
            'ISG SPARROWS POINT (BETHLEHEM STEEL CORP)';...
            'COX CREEK'}
            {'PATAPSCO';...
             'W R GRACE';...
             'ERACHEM'}};

location_name = {'Lower James River-HRSD'
                 'Potomac River-Blue Plains'
                 'Back River'
                 'Patapsco River'};

station_id = 1:size(location,1);
for i = 1:size(location,1)
    tmp = get_station_ij(location{i});
    station_ij(i,:) = tmp(1,:);
end   

figure;
contourf(lon,lat,mask);
hold on;
for i=1:length(station_ij)
    scatter(lon(station_ij(i,1),station_ij(i,2)),lat(station_ij(i,1),station_ij(i,2)),40,'green','filled');
    hold on;
    text(lon(station_ij(i,1),station_ij(i,2)),lat(station_ij(i,1),station_ij(i,2)),location_name{i},'color','r','fontweight','bold');
end

% Fraction in each sigma layer
ZFRACPS = [0.5 0.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

%UNIT OF TBREAK/ps_t: SECS, MINS, HRS, DAYS
TWARPPS = 'DAYS';

FLOW_time =  [datenum(2009,12,1) datenum(2010,1:12,1),datenum(2011,1,1)];
for loc=1:size(location,1)
    %unit of ps: kg/day !!!
    [data_time{loc},~,~,FLOW{loc},~,~,~,~,~,~,~,~,~] = get_ps('2000 Point Sources MD.DE.DC.VA.xlsx',location{loc});
    for i = 1:length(location{loc})
        t = data_time{loc}{i}{3};
        flw = FLOW{loc}{i};
        tu = unique(t);
        for j=1:length(tu)
            tnew(j) = tu(j);
            flw_new(j) = sum(flw(t==tu(j)));
        end
        FLOW_tmp(i,:) = interp1(tnew,flw_new,datenum(2000,1:12,1),'nearest','extrap');
        clear tnew flw_new
    end
    FLOW_mon(loc,:) = sum(FLOW_tmp,1);
    clear FLOW_tmp;
end

FLOW_mean = mean(FLOW_mon,2);
save('WWTP_flow.mat',"station_ij","location_name","location","FLOW_mean");