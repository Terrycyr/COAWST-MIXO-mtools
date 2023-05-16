clear all; close all;

load('non_tidal_bnd_2021.mat');

el1 = s_el_out;
v1 = s_v_out;

load('./Backup/non_tidal_bnd_2021.mat')
el2 = s_el_out;
v2 = s_v_out;

load('./Detided_GOM_HYCOM/non_tidal_bnd_2021.mat')
el3 = s_el_out;
v3 = s_v_out;

el4 = ncread('../../NC_file_generation/WFS_2021_bry.nc','zeta_south')-0.35;
v4 = ncread('../../NC_file_generation/WFS_2021_bry.nc','v_south');

ix = 10;
figure;
subplot(1,2,1)
plot(squeeze(el4(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el2(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el3(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el1(ix,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');

subplot(1,2,2)
plot(squeeze(v4(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v2(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v3(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v1(ix,end,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');

ix = 100;
figure;
subplot(1,2,1)
plot(squeeze(el4(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el2(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el3(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el1(ix,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');

subplot(1,2,2)
plot(squeeze(v4(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v2(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v3(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v1(ix,end,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');

ix = 200;
figure;
subplot(1,2,1)
plot(squeeze(el4(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el2(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el3(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el1(ix,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');

subplot(1,2,2)
plot(squeeze(v4(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v2(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v3(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v1(ix,end,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');

ix = 300;
figure;
subplot(1,2,1)
plot(squeeze(el4(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el2(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el3(ix,:)),'LineWidth',2);
hold on;
plot(squeeze(el1(ix,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');

subplot(1,2,2)
plot(squeeze(v4(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v2(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v3(ix,end,:)),'LineWidth',2);
hold on;
plot(squeeze(v1(ix,end,:)),'LineWidth',2);
hold on;
legend('BND File','Tide', 'Detide','Now');