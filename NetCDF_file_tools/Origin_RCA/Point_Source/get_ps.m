function [data_time,BOD5,DO,FLOW,NH3,NO23,PO4,TKN,TN,TON,TOP,TP,TSS] = get_ps(fn,location)
%get_ps Summary of this function goes here
%   Detailed explanation goes here  

    [dat,~,raw] = xlsread(fn,1);
    st_name = raw(:,1); 

    for i=1:length(location)
        pos = find(strcmpi(st_name,location{i}));
        time = raw(pos,3);
        for j=1:length(time)
            time2(j) = datenum(time{j},'mm/dd/yyyy');
        end

        
        raw2 = raw(pos,:);
        var = raw2(:,4);

        % BOD5
        pos_BOD5 = strcmpi(var,'BOD5');
        data_time{i}{1} = time2(pos_BOD5);
        BOD5{i} = cell2mat(raw2(pos_BOD5,5));

        % DO
        pos_DO = strcmpi(var,'DO');
        data_time{i}{2} = time2(pos_DO);
        DO{i} = cell2mat(raw2(pos_DO,5));

        % FLOW
        pos_FLOW = strcmpi(var,'FLOW');
        data_time{i}{3} = time2(pos_FLOW);
        FLOW{i} = cell2mat(raw2(pos_FLOW,5));

        % NH3
        pos_NH3 = strcmpi(var,'NH3');
        data_time{i}{4} = time2(pos_NH3);
        NH3{i} = cell2mat(raw2(pos_NH3,5));

        % NO23
        pos_NO23 = strcmpi(var,'NO23');
        data_time{i}{5} = time2(pos_NO23);
        NO23{i} = cell2mat(raw2(pos_NO23,5));

        % PO4
        pos_PO4 = strcmpi(var,'PO4');
        data_time{i}{6} = time2(pos_PO4);
        PO4{i} = cell2mat(raw2(pos_PO4,5));

        % TKN
        pos_TKN = strcmpi(var,'TKN');
        data_time{i}{7} = time2(pos_TKN);
        TKN{i} = cell2mat(raw2(pos_TKN,5));

        % TN
        pos_TN = strcmpi(var,'TN');
        data_time{i}{8} = time2(pos_TN);
        TN{i} = cell2mat(raw2(pos_TN,5));

        % TON
        pos_TON = strcmpi(var,'TON');
        data_time{i}{9} = time2(pos_TON);
        TON{i} = cell2mat(raw2(pos_TON,5));

        % TOP
        pos_TOP = strcmpi(var,'TOP');
        data_time{i}{10} = time2(pos_TOP);
        TOP{i} = cell2mat(raw2(pos_TOP,5));

        % TP
        pos_TP = strcmpi(var,'TP');
        data_time{i}{11} = time2(pos_TP);
        TP{i} = cell2mat(raw2(pos_TP,5));

        % TSS
        pos_TSS = strcmpi(var,'TSS');
        data_time{i}{12} = time2(pos_TSS);
        TSS{i} = cell2mat(raw2(pos_TSS,5));
    end
end