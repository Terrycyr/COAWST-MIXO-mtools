clear all;close all;
addpath(path,'C:\Users\cheny\Desktop\EcoHAB\self_functions');

dat = [7 0;
       8.5 0.65;
       11.5 1.1;
       17.5 1.2;
       24 1.18;
       30 0.9];

dat2 = [11 0;
        13 0.24;
        16 1.7;
        19 1.6;
        22.2 2.2;
        25 3.0;
        27.8 2.35;
        31 2.4;
        33.4 1.5;
        36 0];

dat3 = [0., 0.4847206911579969;
    4.909375240341475, 0.9643781885302634;
    10.0488373881611, 1.6479093496039177;
    15.939767221267978, 3.1342451354885026;
    22.024252057323043, 3.8550772938190585];
dat3(:,2) = div2grate(dat3(:,2));

dat4 = [9.876511155838415, 0.37575799190616976;
14.88700395250132, 0.8542295197346119;
19.961183089645818, 1.23679817352382;
24.89702468060271, 1.3111747211576836;
29.896252204774758, 1.275949823403442;
34.93438246994855, 0.014689915328066405];

t = 5:35;
G = 0.275*exp(0.0639*(t));


TOPT1 = 24;
K1C = 1.8;
K1BETA1 = 0.01;
K1BETA2 = 0.012;
temp1 = t;

for i=1:length(temp1)
    if(temp1(i)<=TOPT1)
        GITMAX1(i) = K1C*exp(-K1BETA1*(temp1(i)-TOPT1)^2.);
    elseif(temp1(i)>TOPT1)
        GITMAX1(i) = K1C*exp(-K1BETA2*(TOPT1-temp1(i))^2.);
    end
end

figure('unit','pixel','position',[180   187   1020   284]);
plot(t,GITMAX1,t,G);
hold on;
scatter(dat(:,1),dat(:,2),40,'r','filled');
hold on;
scatter(dat2(:,1),dat2(:,2),40,'g','filled');
hold on;
scatter(dat3(:,1),dat3(:,2),40,'b','filled');
hold on;
scatter(dat4(:,1),dat4(:,2),40,'m','filled');
legend('string',{'Model','Walsh et al.','Rhizosolenia fragilissima'...
                                       ,'Skeletonema Tropicum'...
                                       ,'Skeletonema Costatum (1979)'...
                                       ,'Skeletonema Costatum (2011)'},'location','northwestoutside');
set(gca,'fontsize',15);
set(gcf,'color','w');
xlabel('Temp.');
ylabel('Growth Rate');


dat3 = load('tvsgrowth.txt');
%t = 12:28;
t = 5:35;
G = 0.1*exp(0.0693*(t));

TOPT2 = 27;
K2C = 0.8;
K1BETA1 = 0.01;
K1BETA2 = 0.03;
temp1 = t;

for i=1:length(temp1)
    if(temp1(i)<=TOPT2)
        GITMAX2(i) = K2C*exp(-K1BETA1*(temp1(i)-TOPT2)^2.);
    elseif(temp1(i)>TOPT2)
        GITMAX2(i) = K2C*exp(-K1BETA2*(TOPT2-temp1(i))^2.);
    end
end

figure('unit','pixel','position',[180   187   1020   284]);
plot(t,GITMAX2,t,G);
hold on;
scatter(dat3(:,1),dat3(:,2)*K2C,40,'r','filled');
legend('string',{'Model','Walsh et al.','Literature'},'location','northwestoutside');
set(gca,'fontsize',15);
set(gcf,'color','w');
xlabel('Temp.');
ylabel('Growth Rate');

%dat = load('tvsgrowth.txt');
%t = 12:28;
t = 5:35;
%G = 0.225*exp(0.0693*(t));
G = 1.8*exp(0.0633*(t-27));

TOPT3 = 28;
K3C = 0.9;
K1BETA1 = 0.01;
K1BETA2 = 0.015;
temp1 = t;

for i=1:length(temp1)
    if(temp1(i)<=TOPT3)
        GITMAX3(i) = K3C*exp(-K1BETA1*(temp1(i)-TOPT3)^2.);
    elseif(temp1(i)>TOPT3)
        GITMAX3(i) = K3C*exp(-K1BETA2*(TOPT3-temp1(i))^2.);
    end
end


figure('unit','pixel','position',[180   187   1020   284]);
plot(t,GITMAX3,t,G);
%hold on;
%scatter(dat(:,1),dat(:,2),40,'r','filled');
legend('string',{'Model','Walsh et al.'},'location','northwestoutside');
set(gca,'fontsize',15);
set(gcf,'color','w');
xlabel('Temp.');
ylabel('Growth Rate');


figure('unit','pixel','position',[180   187   1020   284]);
plot(t,GITMAX1,t,GITMAX2,t,GITMAX3);
hold on;
legend('string',{'Diatom','K. brevis','Synechococcus sp.'},'location','northwestoutside');
set(gca,'fontsize',15);
set(gcf,'color','w');
xlabel('Temp.');
ylabel('Growth Rate');

figure('unit','pixel','position',[180   187   1020   284]);
scatter(dat(:,1),dat(:,2),70,'r^','filled');
hold on;
scatter(dat2(:,1),dat2(:,2),70,'gs','filled');
hold on;
scatter(dat3(:,1),dat3(:,2)*K2C,70,'b','filled');
alpha(0.6);
legend('string',{'\it{Rhizosolenia fragilissima}','\it{Skeletonema Tropicum}','\it{K. brevis}'},'location','northwestoutside');
set(gca,'fontsize',15);
set(gcf,'color','w');
xlabel('Temperature');
ylabel('Growth Rate');