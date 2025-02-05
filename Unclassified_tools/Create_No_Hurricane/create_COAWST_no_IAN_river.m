clear all; close all;

river_fname_orig = 'WFS_2022_river.nc';
river_fname_out = 'WFS_2022_river_no_IAN.nc';
hurricane_lanfall = datenum(2022,9,28,20,35,00);
year = 2022;
N = 21;
replace_len_before = 2;
replace_len_after = 15;

var_replace = {'river_transport','river_temp','river_salt','river_mud_01','river_mud_02'};

copyfile(['../',river_fname_orig],['./',river_fname_out]);

river_time = ncread(river_fname_out,'river_time')+datenum(year,1,1);
rep_pos = find((river_time>=hurricane_lanfall-replace_len_before)...
    &(river_time<=hurricane_lanfall+replace_len_after));
river_time2 = river_time;
river_time2(rep_pos) = [];

for ii=1:length(var_replace)
    river_var_orig = ncread(river_fname_out,var_replace{ii});
    river_var2 = river_var_orig;

    if(size(river_var_orig,2) == N)
        river_var2(:,:,rep_pos) = [];
        for i = 1:size(river_var_orig,1)
            for j=1:N
                river_var(i,j,:) = interp1(river_time2,squeeze(river_var2(i,j,:)),river_time);
            end
        end

        var_orig_plt = squeeze(river_var_orig(:,end,:));
        var_plt = squeeze(river_var(:,end,:));
    else
        river_var2(:,rep_pos) = [];
        for i = 1:size(river_var_orig,1)
            river_var(i,:) = interp1(river_time2,river_var2(i,:),river_time);
        end

        var_orig_plt = river_var_orig;
        var_plt = river_var;
    end    

    figure('Units','pixels','Position',[100 100 1300 500])
    ht = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
    for i = 1:size(var_orig_plt,1)
        nexttile;
        plot(river_time,var_plt(i,:),'LineWidth',3);
        hold on
        plot(river_time,var_orig_plt(i,:),'LineWidth',1.5);
        xlim([hurricane_lanfall-30 hurricane_lanfall+30]);
        set(gca,'xtick',[hurricane_lanfall-30:10:hurricane_lanfall+30],...
            'xticklabel',datestr([hurricane_lanfall-30:10:hurricane_lanfall+30],'mm/dd'));
    end
    set(gcf,'color','w');
    title(ht,var_replace{ii}, 'interpreter', 'none');

    ncwrite(river_fname_out,var_replace{ii},river_var);

    clear river_var
end