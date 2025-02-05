clear all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

fn = {'Myakka.xlsx','Peace.xlsx','Caloosahatchee.xlsx'};
st_name = {'3499','PR12','CES06'};

k=0;

%Nutrients

for st = 1:length(st_name)
    [dat,txt,raw] = xlsread(fn{st},1);
    st_all = raw(:,4);
    pos = find(contains(st_all,st_name{st}));
    st_lat(st) = raw{pos(1),7};
    st_lon(st) = raw{pos(1),8};
    for i=1:length(pos)
        date_tmp(i) = xlstime2date(raw{pos(i),10});
        dep_tmp(i) = raw{pos(i),11};
        varcode_tmp{i} = raw{pos(i),13};
        var_QA_tmp{i} = raw{pos(i),18};
        if(strcmp(var_QA_tmp{i},'U'))
            var_tmp(i) = 0.;
        elseif(strcmp(var_QA_tmp{i},'YI'))
            var_tmp(i) = NaN;
        else
            var_tmp(i) = raw{pos(i),16};
        end
    end

%TN
    pos1 = find(contains(varcode_tmp,'TN_mgl'));
    pos2 = find(contains(varcode_tmp,'TN_ugl'));

    tmp1 = date_tmp(pos1);
    tmp2 = date_tmp(pos2);
    st_tn{st,1} = [tmp1,tmp2];

    tmp1 = dep_tmp(pos1);
    tmp2 = dep_tmp(pos2);
    st_tn{st,2} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2)/1000;
    st_tn{st,3} = [tmp1,tmp2];

    tmp1 = varcode_tmp(pos1);
    tmp2 = varcode_tmp(pos2);
    st_tn{st,4} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    st_tn{st,5} = [tmp1,tmp2];
    
%Ammonium
    pos1 = find(contains(varcode_tmp,'NH4_mgl'));
    pos2 = find(contains(varcode_tmp,'NH4_ugl'));
    pos3 = find(contains(varcode_tmp,'NH3_N_mgl'));
    pos4 = find(contains(varcode_tmp,'NH3_N_ugl'));

    tmp1 = date_tmp(pos1);
    tmp2 = date_tmp(pos2);
    tmp3 = date_tmp(pos3);
    tmp4 = date_tmp(pos4);
    st_nh4{st,1} = [tmp1,tmp2,tmp3,tmp4];

    tmp1 = dep_tmp(pos1);
    tmp2 = dep_tmp(pos2);
    tmp3 = dep_tmp(pos3);
    tmp4 = dep_tmp(pos4);
    st_nh4{st,2} = [tmp1,tmp2,tmp3,tmp4];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2)/1000;
    tmp3 = var_tmp(pos3);
    tmp4 = var_tmp(pos4)/1000;
    st_nh4{st,3} = [tmp1,tmp2,tmp3,tmp4];

    tmp1 = varcode_tmp(pos1);
    tmp2 = varcode_tmp(pos2);
    tmp3 = varcode_tmp(pos3);
    tmp4 = varcode_tmp(pos4);
    st_nh4{st,4} = [tmp1,tmp2,tmp3,tmp4];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    tmp3 = var_tmp(pos3);
    tmp4 = var_tmp(pos4);
    st_nh4{st,5} = [tmp1,tmp2,tmp3,tmp4];

%Nitrate
    pos1 = find(contains(varcode_tmp,'NOx_mgl'));
    pos2 = find(contains(varcode_tmp,'NOx_ugl'));

    tmp1 = date_tmp(pos1);
    tmp2 = date_tmp(pos2);
    st_no23{st,1} = [tmp1,tmp2];

    tmp1 = dep_tmp(pos1);
    tmp2 = dep_tmp(pos2);
    st_no23{st,2} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2)/1000;
    st_no23{st,3} = [tmp1,tmp2];

    tmp1 = varcode_tmp(pos1);
    tmp2 = varcode_tmp(pos2);
    st_no23{st,4} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    st_no23{st,5} = [tmp1,tmp2];

%ON
    pos1 = find(contains(varcode_tmp,'norg_ugl'));

    tmp1 = date_tmp(pos1);
    st_on{st,1} = [tmp1];

    tmp1 = dep_tmp(pos1);
    st_on{st,2} = [tmp1];

    tmp1 = var_tmp(pos1);
    st_on{st,3} = [tmp1/1000];

    tmp1 = varcode_tmp(pos1);
    st_on{st,4} = [tmp1];

    tmp1 = var_tmp(pos1);
    st_on{st,5} = [tmp1];

%ON-Kjeldahl
    if(isempty(st_on{st,3})&&~isempty(st_nh4{st,3}))
        pos1 = find(contains(varcode_tmp,'TKN_mgl'));
        pos2 = find(contains(varcode_tmp,'TKN_ugl'));
        tmp1 = date_tmp(pos1);
        tmp2 = date_tmp(pos2);
        tkn{1} = [tmp1,tmp2];

        tmp1 = dep_tmp(pos1);
        tmp2 = dep_tmp(pos2);
        tkn{2} = [tmp1,tmp2];

        tmp1 = var_tmp(pos1);
        tmp2 = var_tmp(pos2)/1000;
        tkn{3} = [tmp1,tmp2];

        tmp1 = varcode_tmp(pos1);
        tmp2 = varcode_tmp(pos2);
        tkn{4} = [tmp1,tmp2];

        tmp1 = var_tmp(pos1);
        tmp2 = var_tmp(pos2);
        tkn{5} = [tmp1,tmp2];

        if(~isempty(tkn{3}))
            x1 = st_nh4{st,1};
            y1 = st_nh4{st,3};
            x2 = tkn{1};
            y2 = tkn{3};
            [~,~,tcate] = unique(x1);
            for t=1:max(tcate)
                x12(t) = mean(x1(tcate==t));
                y12(t) = mean(y1(tcate==t));
            end

            on_date = tkn{1};
            on_dep = tkn{2};

            if(length(x12)>1)
                on_tmp = max(y2-interp1(x12,y12,x2),0);

                if(~isempty(on_tmp(~isnan(on_tmp))))
                    st_on{st,1} = on_date(~isnan(on_tmp));
                    st_on{st,2} = on_dep(~isnan(on_tmp));
                    st_on{st,3} = on_tmp(~isnan(on_tmp));
                    st_on{st,4} = repmat({'on_mgl'},1,length(st_on{st,3}));
                    st_on{st,5} = 'Calculated from Kjeldahl';
                end
            end
        end
    end

    clear x12 y12

%TP
    pos1 = find(contains(varcode_tmp,'TP_mgl'));
    pos2 = find(contains(varcode_tmp,'TP_ugl'));

    tmp1 = date_tmp(pos1);
    tmp2 = date_tmp(pos2);
    st_tp{st,1} = [tmp1,tmp2];

    tmp1 = dep_tmp(pos1);
    tmp2 = dep_tmp(pos2);
    st_tp{st,2} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2)/1000;
    st_tp{st,3} = [tmp1,tmp2];

    tmp1 = varcode_tmp(pos1);
    tmp2 = varcode_tmp(pos2);
    st_tp{st,4} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    st_tp{st,5} = [tmp1,tmp2];

%Phospahte
    pos1 = find(contains(varcode_tmp,'OP_mgl'));
    pos2 = find(contains(varcode_tmp,'OP_ugl'));

    tmp1 = date_tmp(pos1);
    tmp2 = date_tmp(pos2);
    st_po4{st,1} = [tmp1,tmp2];

    tmp1 = dep_tmp(pos1);
    tmp2 = dep_tmp(pos2);
    st_po4{st,2} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2)/1000;
    st_po4{st,3} = [tmp1,tmp2];

    tmp1 = varcode_tmp(pos1);
    tmp2 = varcode_tmp(pos2);
    st_po4{st,4} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    st_po4{st,5} = [tmp1,tmp2];

%Silica
    pos1 = find(contains(varcode_tmp,'Si_mgl'));
    pos2 = find(contains(varcode_tmp,'Si_ugl'));
    pos3 = find(contains(varcode_tmp,'SiO2_ugl'));

    tmp1 = date_tmp(pos1);
    tmp2 = date_tmp(pos2);
    tmp3 = date_tmp(pos3);
    st_si{st,1} = [tmp1,tmp2,tmp3];

    tmp1 = dep_tmp(pos1);
    tmp2 = dep_tmp(pos2);
    tmp3 = dep_tmp(pos3);
    st_si{st,2} = [tmp1,tmp2,tmp3];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2)/1000;
    tmp3 = var_tmp(pos3)/1000/60*28;
    st_si{st,3} = [tmp1,tmp2,tmp3];

    tmp1 = varcode_tmp(pos1);
    tmp2 = varcode_tmp(pos2);
    tmp3 = varcode_tmp(pos3);
    st_si{st,4} = [tmp1,tmp2,tmp3];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    tmp3 = var_tmp(pos3);
    st_si{st,5} = [tmp1,tmp2,tmp3];

%Chla
    pos1 = find(contains(varcode_tmp,'Chla_ugl'));
    pos2 = find(contains(varcode_tmp,'ChlaC_ugl'));

    tmp1 = date_tmp(pos1);
    tmp2 = date_tmp(pos2);
    st_chla{st,1} = [tmp1,tmp2];

    tmp1 = dep_tmp(pos1);
    tmp2 = dep_tmp(pos2);
    st_chla{st,2} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    st_chla{st,3} = [tmp1,tmp2];

    tmp1 = varcode_tmp(pos1);
    tmp2 = varcode_tmp(pos2);
    st_chla{st,4} = [tmp1,tmp2];

    tmp1 = var_tmp(pos1);
    tmp2 = var_tmp(pos2);
    st_chla{st,5} = [tmp1,tmp2];

    clear date_tmp dep_tmp varcode_tmp var_tmp
end

%River discharge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_dir = 'C:\Users\cheny\Desktop\EcoHAB\USGS_River\';
r(1) = load([r_dir,'USGS_02298880_MYAKKA.mat']);
r(2) = load([r_dir,'USGS_02296750_PEACE.mat']);
r(3) = load([r_dir,'USGS_02292900_CALOOSAHATCHEE.mat']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
river_flow_date = datenum(2000,1,1):datenum(2024,1,1);
river_flow = get_nwis_timeseries(r,river_flow_date,'Mean river discharge');


save('2000-2024-3River-Data.mat',"st_name","st_lon","st_lat","st_tn","st_no23"...
    ,"st_nh4","st_on","st_tp","st_po4","st_si","st_chla","river_flow_date"...
    ,"river_flow");