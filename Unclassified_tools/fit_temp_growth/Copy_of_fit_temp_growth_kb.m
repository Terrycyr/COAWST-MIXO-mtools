clear all; close all

t = 10:35;
% 

TOPT1 = 24;
K1C = 2.4;
K1BETA1 = 0.008;
K1BETA2 = 0.01;
for i=1:length(t)
    if(t(i)<=TOPT1)
        GITMAX1(i) = K1C*exp(-K1BETA1*(t(i)-TOPT1)^2.);
    else
        GITMAX1(i) = K1C*exp(-K1BETA2*(TOPT1-t(i))^2.);
    end
end

TOPT12 = 24;
K1C2 = 1.8;
K1BETA12 = 0.008;
K1BETA22 = 0.01;
for i=1:length(t)
    if(t(i)<=TOPT1)
        GITMAX12(i) = K1C2*exp(-K1BETA12*(t(i)-TOPT12)^2.);
    else
        GITMAX12(i) = K1C2*exp(-K1BETA22*(TOPT12-t(i))^2.);
    end
end

figure('unit','pixel','position',[180   187   1020   284]);
plot(t,GITMAX1);
hold on;
plot(t,GITMAX12);
hold on;
legend('2021','2022');
set(gca,'fontsize',15);
set(gcf,'color','w');
xlabel('Temp.');
ylabel('Growth Rate');
