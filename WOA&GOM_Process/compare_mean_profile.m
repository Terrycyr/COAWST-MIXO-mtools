clear all;
load('./ocean_nutrients_mean_profile2.mat');
no3_m2 = no3_m;
po4_m2 = po4_m;
load('./ocean_nutrients_mean_profile.mat');

NEGOM_PO4 = xlsread('./NEGOM_nutrients.xlsx',1);
NEGOM_NO3 = xlsread('./NEGOM_nutrients.xlsx',2);

figure;
plot(NEGOM_NO3(:,1)/1.025,NEGOM_NO3(:,2)*-1,no3_m, woa_dep*-1,no3_m2, woa_dep*-1,'LineWidth',1.5);set(gcf,'color','w');
legend('NEGOM-2000-04','WOA18','WOA18_modified');
ylim([-500 0]);
grid on;
xlabel('\muM NO3');
ylabel('Depth m');

figure;
plot(NEGOM_NO3(:,1)/1.025,NEGOM_NO3(:,2)*-1,no3_m, woa_dep*-1,no3_m2, woa_dep*-1,'LineWidth',1.5);
set(gcf,'color','w');
legend('NEGOM-2000-04','WOA18','WOA18_modified');
ylim([-1000 0]);
grid on;
xlabel('\muM NO3');
ylabel('Depth m');

figure;
plot(NEGOM_PO4(:,1)/1.025,NEGOM_PO4(:,2)*-1,po4_m, woa_dep*-1,po4_m, woa_dep*-1,'LineWidth',1.5);
set(gcf,'color','w');
legend('NEGOM-2000-04','WOA18','WOA18_modified');
ylim([-500 0]);
grid on;
xlabel('\muM po4');
ylabel('Depth m');

figure;
plot(NEGOM_PO4(:,1)/1.025,NEGOM_PO4(:,2)*-1,po4_m, woa_dep*-1,po4_m, woa_dep*-1,'LineWidth',1.5);
set(gcf,'color','w');
legend('NEGOM-2000-04','WOA18','WOA18_modified');
ylim([-1000 0]);
grid on;
xlabel('\muM po4');
ylabel('Depth m');