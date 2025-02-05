clear all; close all;

year = 2022;
hurricane_lanfall = datenum(2022,9,28,20,35,00)-datenum(year,1,1);
N = 21;
replace_len_before = 3;
replace_len_after = 5;


nontidal_fname_orig = ['non_tidal_bnd_',num2str(year),'.mat'];
nontidal_fname_out = ['non_tidal_bnd_',num2str(year),'_no_IAN.mat'];

copyfile(['../../Non_tidal_component_preprocessing/HYCOM/',nontidal_fname_orig],['./',nontidal_fname_out]);

var_list = who('-file',nontidal_fname_out);

load(nontidal_fname_out,'date_out2d');
load(nontidal_fname_out,'date_out3d');

nontidal_time_2d = date_out2d-datenum(year,1,1);
nontidal_time_3d = date_out3d-datenum(year,1,1);

rep_pos_2d = find((nontidal_time_2d>=hurricane_lanfall-replace_len_before)...
    &(nontidal_time_2d<=hurricane_lanfall+replace_len_after));
rep_pos_3d = find((nontidal_time_3d>=hurricane_lanfall-replace_len_before)...
    &(nontidal_time_3d<=hurricane_lanfall+replace_len_after));

rep_pos_2d2 = rep_pos_2d;
rep_pos_2d2(1:6) = [];
rep_pos_2d2(end-5:end) = [];
rep_pos_3d2 = rep_pos_3d;
rep_pos_3d2(1) = [];
rep_pos_3d2(end) = [];

nontidal_time_2d2 = nontidal_time_2d;
nontidal_time_2d2(rep_pos_2d) = [];

nontidal_time_3d2 = nontidal_time_3d;
nontidal_time_3d2(rep_pos_3d) = [];

for i = 1:size(var_list,1)
    if(isempty(strfind(var_list{i},'date')))
        load(nontidal_fname_out,var_list{i});
        eval(['tmp_orig=',var_list{i},';']);
        ndim = sum(size(tmp_orig)>1);

        if(ndim==2)
            [a,b] = size(tmp_orig);
            
            for j = 1:a
                tmp2 = smooth(tmp_orig(j,:),24);
                tmp2(rep_pos_2d) = [];
                tmp(j,:) = tmp_orig(j,:);
                tmp(j,rep_pos_2d2) = interp1(nontidal_time_2d2,tmp2,nontidal_time_2d(rep_pos_2d2));
            end

            var_plt = tmp;
            var_orig_plt = tmp_orig;
            nontidal_time = nontidal_time_2d;
        elseif(ndim==3)
            [a,b,c] = size(tmp_orig);

            for j = 1:a
                for k=1:b
                    tmp2 = smooth(squeeze(tmp_orig(j,k,:)),4);
                    tmp2(rep_pos_3d) = [];
                    tmp(j,k,:) = tmp_orig(j,k,:);
                    tmp(j,k,rep_pos_3d2) = interp1(nontidal_time_3d2,tmp2,nontidal_time_3d(rep_pos_3d2));
                end
            end

            var_plt = squeeze(tmp(:,end,:));
            var_orig_plt = squeeze(tmp_orig(:,end,:));
            nontidal_time = nontidal_time_3d;
        end

        eval([var_list{i},'=tmp;']);

        figure('Units','pixels','Position',[100 100 1300 500])
        ht = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
        plt_pos = round(linspace(1,a,9));
        for j = 1:9
            nexttile;
            plot(nontidal_time,var_plt(plt_pos(j),:),'LineWidth',3);
            hold on
            plot(nontidal_time,var_orig_plt(plt_pos(j),:),'LineWidth',1.5);
            xlim([hurricane_lanfall-30 hurricane_lanfall+30]);
            set(gca,'xtick',[hurricane_lanfall-30:10:hurricane_lanfall+30],...
                'xticklabel',datestr([hurricane_lanfall-30:10:hurricane_lanfall+30]+datenum(year,1,1),'mm/dd'));
        end
        set(gcf,'color','w');
        title(ht,var_list{i}, 'interpreter', 'none');

        %eval(['save(nontidal_fname_out,',var_list{i},',''-append'');']);
        clear tmp tmp2
        eval(['clear ',var_list{i},';']);
    end
end