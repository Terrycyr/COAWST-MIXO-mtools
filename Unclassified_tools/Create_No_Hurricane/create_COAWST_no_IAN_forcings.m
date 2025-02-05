clear all; close all;

year = 2022;
hurricane_lanfall = datenum(2022,9,28,20,35,00)-datenum(year,1,1);
N = 21;
replace_len_before = 2;
replace_len_after = 3.5;

var_fname = {'Pair','Uwind','Vwind'};
var_replace = {'Pair','Uwind','Vwind'};
methods =   [    1       2       2];

plt_i = [10 100 200 300 10 100 200 300  10 100 200 300 350];
plt_j = [10  10  10  10 90  90  90  90 180 180 180 180 102];

scale = 0.1;

for nsec = 1
    for ii=1:length(var_replace)

        forcing_fname_orig = ['WFS_2022_',var_fname{ii},num2str(nsec),'.nc'];
        forcing_fname_out = ['WFS_2022_',var_fname{ii},num2str(nsec),'_no_IAN.nc'];

        copyfile(['../',forcing_fname_orig],['./',forcing_fname_out]);

        forcing_time = ncread(forcing_fname_out,'time');
        rep_pos = find((forcing_time>=hurricane_lanfall-replace_len_before)...
            &(forcing_time<=hurricane_lanfall+replace_len_after));

        forcing_time2 = forcing_time;
        forcing_time2(rep_pos) = [];


        forcing_var_orig = ncread(forcing_fname_out,var_replace{ii});

        if(methods(ii)==1)
            forcing_var2 = forcing_var_orig;
            forcing_var2(:,:,rep_pos) = [];
            for i = 1:size(forcing_var_orig,1)
                for j=1:size(forcing_var_orig,2)
                    tmp = squeeze(forcing_var2(i,j,:));
                    tmp = smooth(tmp,4);
                    tmp2 =  interp1(forcing_time2,tmp,forcing_time);
                    forcing_var(i,j,:) = forcing_var_orig(i,j,:);
                    forcing_var(i,j,rep_pos) = tmp2(rep_pos);
                end
            end
        elseif(methods(ii)==2)
            forcing_var = forcing_var_orig;
            forcing_var(:,:,rep_pos) = forcing_var_orig(:,:,rep_pos)*scale;
        end


        for i=1:length(plt_i)
            var_orig_plt(i,:) = squeeze(forcing_var_orig(plt_i(i),plt_j(i),:));
            var_plt(i,:) = squeeze(forcing_var(plt_i(i),plt_j(i),:));
        end


        figure('Units','pixels','Position',[100 100 1300 500])
        ht = tiledlayout(3,5,'TileSpacing','compact','Padding','compact');
        for i = 1:length(plt_i)
            nexttile;
            plot(forcing_time,var_plt(i,:),'LineWidth',3);
            hold on
            plot(forcing_time,var_orig_plt(i,:),'LineWidth',1.5);
            xlim([hurricane_lanfall-30 hurricane_lanfall+30]);
            set(gca,'xtick',[hurricane_lanfall-30:10:hurricane_lanfall+30],...
                'xticklabel',datestr([hurricane_lanfall-30:10:hurricane_lanfall+30]+datenum(year,1,1),'mm/dd'));
        end
        set(gcf,'color','w');
        title(ht,var_replace{ii}, 'interpreter', 'none');

        ncwrite(forcing_fname_out,var_replace{ii},forcing_var);

        clear forcing_var forcing_time
    end
end