function [calflow,caldat] = create_loadest_file(r, rname, varname, var_flag, n_river, out_date)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
year = datevec(out_date(1));
year = year(1);
load(strcat('WA_',num2str(year),'.mat'));
out_date2 = datenum(year,1,1,0,0,0):datenum(year,12,31,24,0,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calflow = cell(n_river,length(varname));
caldat = cell(n_river,length(varname));
for river = 1:n_river

    for i=1:length(varname)
        if(var_flag(river,i)==1)
            out_path = ['.\',rname{river},'\',varname{i},'\'];
            if(~exist(out_path))
                mkdir(out_path);
            end

            copyfile('./loadest.exe',out_path);

            out_date0 = nutri_dat{river,i}(:,1);
            river_flow = get_nwis_timeseries(r(river),out_date0,'Mean river discharge');
            calib_f = river_flow*35.3147;
            calib_v = nutri_dat{river,i}(:,2);
            calib_v(calib_f<=0) = [];
            calib_f(calib_f<=0) = [];
            calib_date = datestr(out_date0,'yyyymmdd');
            calib_time = datestr(out_date0,'    HHMM');

            calflow{river,i} = calib_f;
            caldat{river,i} = calib_v;


            river_flow_e = get_nwis_timeseries(r(river),out_date2,'Mean river discharge');
            est_f = max(1,river_flow_e*35.3147);
            est_date = datestr(out_date2,'yyyymmdd');
            est_time = datestr(out_date2,'    HHMM');


            fid = fopen([out_path,'header.inp'],'wt+');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#  LOADEST Header File');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n',[rname{river} ' River at Florida']);
            fprintf(fid,'%s\n','1             |     PRTOPT (col.1-5)');
            fprintf(fid,'%s\n','3             |     SEOPT (col.1-5)');
            fprintf(fid,'%s\n','0             |     LDOPT (col. 1-5)');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','# model number, MODNO (col.1-5)');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','1');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','# number of constituents, NCONST (col.1-5)');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','1');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#  Unit flags and constituent names, for I=1,NCONST');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#                                            Unit Flags');
            fprintf(fid,'%s\n','#CNAME                                       Conc Load');
            fprintf(fid,'%s\n','#                                            |    |');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n',[varname{i},'                                           1    1']);
            fclose(fid);

            fid = fopen([out_path,'control.inp'],'wt+');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#    LOADEST control file');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#    line              name of the:');
            fprintf(fid,'%s\n','#    ----              --------------');
            fprintf(fid,'%s\n','#     1                header file');
            fprintf(fid,'%s\n','#     2                calibration file');
            fprintf(fid,'%s\n','#     3                estimation file');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','header.inp');
            fprintf(fid,'%s\n','calib.inp');
            fprintf(fid,'%s\n','est.inp');
            fclose(fid);

            fid = fopen([out_path,'calib.inp'],'wt+');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#  LOADEST Calibration File');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#CDATE      CTIME     CFLOW     CCONC');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            for j=1:length(calib_f)
                fprintf(fid,'%s%s      %-10.3f%-10.3f\n',calib_date(j,:),calib_time(j,:),calib_f(j),calib_v(j));
            end
            fclose(fid);

            fid = fopen([out_path,'est.inp'],'wt+');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#  LOADEST Estimation File');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','#  Number of observations per day, NOBSPD (col. 1-5)');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','1');
            fprintf(fid,'%s\n','######################################################################');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','# EDATE     ETIME     EFLOW');
            fprintf(fid,'%s\n','#');
            fprintf(fid,'%s\n','######################################################################');
            for j=1:length(river_flow_e)
                fprintf(fid,'%s%s      %-10.0f\n',est_date(j,:),est_time(j,:),est_f(j));
            end
            fclose(fid);

        end
    end
end
end

