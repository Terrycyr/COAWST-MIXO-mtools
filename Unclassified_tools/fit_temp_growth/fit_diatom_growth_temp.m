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
K1C = 2.0;
K1BETA1 = 0.007;
K1BETA2 = 0.012;
temp1 = t;

for i=1:length(temp1)
    if(temp1(i)<=TOPT1)
        GITMAX1(i) = K1C*exp(-K1BETA1*(temp1(i)-TOPT1)^2.);
    elseif(temp1(i)>TOPT1)
        GITMAX1(i) = K1C*exp(-K1BETA2*(TOPT1-temp1(i))^2.);
    end
end


TOPT1 = 24;
K1C2 = 1.8;
K1BETA1 = 0.007;
K1BETA2 = 0.012;
temp1 = t;

for i=1:length(temp1)
    if(temp1(i)<=TOPT1)
        GITMAX12(i) = K1C2*exp(-K1BETA1*(temp1(i)-TOPT1)^2.);
    elseif(temp1(i)>TOPT1)
        GITMAX12(i) = K1C2*exp(-K1BETA2*(TOPT1-temp1(i))^2.);
    end
end

t1 = 25;
g1 = K1C*exp(-K1BETA2*(TOPT1-t1)^2.);
t2 = 32;
g2 = K1C2*exp(-K1BETA2*(TOPT1-t2)^2.);
x = [t1 t2];
y = [g1 g2];


TOPT1 = 24;
K1C3 = 2.0;
K1BETA1 = 0.010;
K1BETA2 = 0.015;
temp1 = t;

for i=1:length(temp1)
    if(temp1(i)<=TOPT1)
        GITMAX13(i) = K1C3*exp(-K1BETA1*(temp1(i)-TOPT1)^2.);
    elseif(temp1(i)>TOPT1)
        GITMAX13(i) = K1C3*exp(-K1BETA2*(TOPT1-temp1(i))^2.);
    end
end

figure('unit','pixel','position',[180   187   1020   284]);
plot(t,GITMAX1);
hold on;
plot(t,GITMAX12);
hold on;
plot(t,GITMAX13);
hold on;
scatter(x,y,'black','filled');
legend('2021','2022','New','2021&2022');
set(gca,'fontsize',15);
set(gcf,'color','w');
xlabel('Temp.');
ylabel('Growth Rate');




