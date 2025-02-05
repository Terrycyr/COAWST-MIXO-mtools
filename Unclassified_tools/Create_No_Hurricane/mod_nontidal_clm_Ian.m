clear all; close all;

year = 2022;
hurricane_lanfall = datenum(2022,9,28,20,35,00)-datenum(year,1,1);
N = 21;
replace_len_before = 3;
replace_len_after = 6;


nontidal_fname_orig = ['non_tidal_clm_',num2str(year),'.mat'];
nontidal_fname_out = ['non_tidal_clm_',num2str(year),'_no_IAN.mat'];

copyfile(['../../Non_tidal_component_preprocessing/HYCOM/',nontidal_fname_orig],['./',nontidal_fname_out]);

var_list = who('-file',nontidal_fname_out);

load(nontidal_fname_out,'date_clm');

nontidal_time_clm = date_clm;

rep_pos_clm = find((nontidal_time_clm>=hurricane_lanfall-replace_len_before)...
    &(nontidal_time_clm<=hurricane_lanfall+replace_len_after));

rep_pos_clm2 = rep_pos_clm;
rep_pos_clm2(1) = [];
rep_pos_clm2(end) = [];

nontidal_time_clm2 = nontidal_time_clm;
nontidal_time_clm2(rep_pos_clm) = [];

plt_i = [10 100 200 300 10 100 200 300  10 100 200 300];
plt_j = [ 1   1   1   1  4   4   4   4   8   8   8   8];

for i = 1:size(var_list,1)
    if(isempty(strfind(var_list{i},'date')))
        load(nontidal_fname_out,var_list{i});
        eval(['tmp_orig=',var_list{i},';']);
        ndim = sum(size(tmp_orig)>1);

        if(ndim==3)
            [a,b,c] = size(tmp_orig);
            
            for j = 1:a
                for k=1:b
                    tmp2 = smooth(tmp_orig(j,k,:),1);
                    tmp2(rep_pos_clm) = [];
                    tmp(j,k,:) = tmp_orig(j,k,:);
                    tmp(j,k,rep_pos_clm2) = interp1(nontidal_time_clm2,tmp2,nontidal_time_clm(rep_pos_clm2));
                end
            end

            for j=1:length(plt_i)
                var_plt(j,:) = squeeze(tmp(plt_i(j),plt_j(j),:));
                var_orig_plt(j,:) = squeeze(tmp_orig(plt_i(j),plt_j(j),:));
                nontidal_time = nontidal_time_clm;
            end

        elseif(ndim==4)

            [a,b,c,d] = size(tmp_orig);

            for j = 1:a
                for k=1:b
                    for l=1:c
                        tmp2 = smooth(squeeze(tmp_orig(j,k,l,:)),4);
                        tmp2(rep_pos_clm) = [];
                        tmp(j,k,l,:) = tmp_orig(j,k,l,:);
                        tmp(j,k,l,rep_pos_clm2) = interp1(nontidal_time_clm2,tmp2,nontidal_time_clm(rep_pos_clm2));
                    end
                end
            end

            for j=1:length(plt_i)
                var_plt(j,:) = squeeze(tmp(plt_i(j),plt_j(j),end,:));
                var_orig_plt(j,:) = squeeze(tmp_orig(plt_i(j),plt_j(j),end,:));
                nontidal_time = nontidal_time_clm;
            end
        end

        eval([var_list{i},'=tmp;']);

        figure('Units','pixels','Position',[100 100 1300 500])
        ht = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');

        for j = 1:length(plt_i)
            nexttile;
            plot(nontidal_time,var_plt(j,:),'LineWidth',3);
            hold on
            plot(nontidal_time,var_orig_plt(j,:),'LineWidth',1.5);
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