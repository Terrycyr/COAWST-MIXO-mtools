function [] = gen_river_nutrients(grd,varname,var_flag,n_river,out_date,fn_WA,fn_loadest,fn_out)

% params.
fn = grd;
N= length(ncread(fn,'Cs_r'));
Vtransform = ncread(fn,'Vtransform');                         
Vstretching = ncread(fn,'Vstretching');
THETA_S = ncread(fn,'theta_s');                     
THETA_B = ncread(fn,'theta_b');                     
TCLINE = ncread(fn,'Tcline');
hc = ncread(fn,'hc');
year = datevec(out_date(1));
year = year(1);
nz = N;

load(fn_loadest);
load(fn_WA);

%DISSOLVED & PARTICULATE
OCDP_r = 0.84;
ONDP_r = 0.60;
OPDP_r = 0.22;

%LABILE & REFRACTORY
OCLR_r = 0.35;
ONLR_r = 0.35;
OPLR_r = 0.35;
      
TIME = out_date-out_date(2);


for river = 1:n_river
    for var = 1:length(varname)
        if(var_flag(river,var)==0)
                try
                    tmp = interp1(nutri_dat{river,var}(:,1),nutri_dat{river,var}(:,2),...
                        out_date);
                    tmp(isnan(tmp)) = interp1(nutri_dat{river,var}(:,1),nutri_dat{river,var}(:,2),...
                        out_date(isnan(tmp)),'nearest','extrap');
                    nutrient1{river}(var,:) = tmp;
                catch
                    nutrient1{river}(var,:) = 0.0;
                end
        else
            nutrient1{river}(var,:) = interp1(datenum(year,1,1,0,0,0):datenum(year,12,31,24,0,0),...
                result_conc{river,var}(:,2),out_date);
        end
    end
end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% nutrient1{1} = nutrient1{1}./10;
% nutrient1{2} = nutrient1{2}./10;
% nutrient1{3} = nutrient1{3}./10;
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%SI
river_DO = 7*ones(n_river,nz,length(TIME));
river_SAL = 0*ones(n_river,nz,length(TIME));
for river = 1:n_river
    OC_ON = (nutrient1{river}(8,:)/12)./(nutrient1{river}(3,:)/14);
    OC_OP = (nutrient1{river}(8,:)/12)./(nutrient1{river}(9,:)/31);
    
    LON = 1-(OC_ON-3.6)/(33.6-3.6);
    LON(LON>1) = 1;
    LON(LON<0) = 0;
    LON_RATIO(river,:) = LON;
    
    LOP = 1-(OC_OP-22.4)/(662.9-22.4);
    LOP(LOP>1) = 1;
    LOP(LOP<0) = 0;
    LOP_RATIO(river,:) = LOP;
    
    for z = 1:nz
        river_LDOC(river,z,:) = OCDP_r*OCLR_r*nutrient1{river}(8,:); 
        river_LDON(river,z,:) = ONDP_r*ONLR_r*nutrient1{river}(3,:);
        river_LDOP(river,z,:) = OPDP_r*OPLR_r*nutrient1{river}(9,:);
        river_LPOC(river,z,:) = (1-OCDP_r)*OCLR_r*nutrient1{river}(8,:);
        river_LPON(river,z,:) = (1-ONDP_r)*ONLR_r*nutrient1{river}(3,:);
        river_LPOP(river,z,:) = (1-OPDP_r)*OPLR_r*nutrient1{river}(9,:);
        river_NH4T(river,z,:) = nutrient1{river}(2,:);
        river_NO23(river,z,:) = nutrient1{river}(1,:);
        river_PO4T(river,z,:) = nutrient1{river}(5,:);
        river_RDOC(river,z,:) = OCDP_r*(1-OCLR_r)*nutrient1{river}(8,:); 
        river_RDON(river,z,:) = ONDP_r*(1-ONLR_r)*nutrient1{river}(3,:);
        river_RDOP(river,z,:) = OPDP_r*(1-OPLR_r)*nutrient1{river}(9,:);
        river_RPOC(river,z,:) = (1-OCDP_r)*(1-OCLR_r)*nutrient1{river}(8,:);
        river_RPON(river,z,:) = (1-OCDP_r)*(1-OCLR_r)*nutrient1{river}(3,:);
        river_RPOP(river,z,:) = (1-OCDP_r)*(1-OCLR_r)*nutrient1{river}(9,:);
        river_SIT(river,z,:) = nutrient1{river}(7,:);     
        if(sum(contains(varname,'CHLA'))>0)
            river_CHLA(river,z,:) = nutrient1{river}(10,:);
        end
    end
end

save(fn_out,'river_LDOC','river_LDON','river_LDOP','river_LPOC','river_LPON','river_LPOP',...
                            'river_NH4T','river_NO23','river_PO4T','river_RDOC','river_RDON','river_RDOP',...
                            'river_RPOC','river_RPON','river_RPOP','river_SIT','river_SAL','river_DO');
if(sum(contains(varname,'CHLA'))>0)
    save(fn_out,'river_CHLA','-append');
end
        
end


