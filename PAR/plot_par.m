clear all; close all;
load('PAR.mat');
time = datenum(2001,1,1)+[0:365];

figure('unit','pixel','position',[180   187   486*1.2   284*1.2]);
plot(time,ITOTSF);
set(gcf,'color','w');
set(gca,'xtick',datenum(2001,1:12,1),'xticklabel',datestr(datenum(2001,1:12,1),'mmm'));
ylabel('PAR (\muE m^-^2 s^-1');
grid on;