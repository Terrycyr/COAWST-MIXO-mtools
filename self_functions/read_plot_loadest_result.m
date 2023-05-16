function [result_flow,result_conc] = read_plot_loadest_result(rname, varname, var_flag, n_river, out_date, fig_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
year = datevec(out_date(1));
year = year(1);
load(strcat('WA_',num2str(year),'.mat'));
load(strcat('calib_',num2str(year),'.mat'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result_flow = cell(n_river,length(varname));
result_conc= cell(n_river,length(varname));
for river=1:n_river
for i=1:length(varname)
    if(var_flag(river,i)==1)
    out_path = ['./',rname{river},'/',varname{i},'/'];
    if(~exist(out_path))
        mkdir(out_path);
    end

    fid = fopen([out_path,varname{i},'.ind'],'rt');
    flag = 0;
    while(flag==0)
        dat = fgetl(fid);
        if(~isempty(dat)&strcmp(dat(10),'-'))
            flag = 1;
        end
    end
    for j=1:366
        result(j,:) = str2num(fgetl(fid));
        result_flow{river,i}(j) = result(j,3);
        tmp = result(j,4:6)*1e6/(result(j,3)*28.316847*24*3600); 
        u_lim = 1.2*max(caldat{river,i});
        tmp(tmp>u_lim) = u_lim;
        result_conc{river,i}(j,:) = tmp;
    end   
    end
end

ax(river,:) = plot_var_flow2(fig_path,calflow(river,:),caldat(river,:)...
    ,result_flow(river,:),result_conc(river,:)...
    ,varname,['Var_dischage_result_',rname{river}],[3 3],'ft^3/s','mg/L',[], year);
end

end
    