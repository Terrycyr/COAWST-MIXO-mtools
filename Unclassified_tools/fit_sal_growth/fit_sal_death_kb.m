clear all; close all;

sal1 = [15 17.5 20 25 30 35];
dat1 = [7 0;
        5 6+10;
        0 2+17;
        0 2+8;
        0 1+18;
        0 2+14];

sal2 = [20 22.5];
dat2 = [6 4;
        0 2+7];

sal3 = [17 18.5];
dat3 = [4 0;
        0 3];

sal4 = [17 18.5];
dat4 = [4 0;
        0 2];

sal = [sal1,sal2,sal3,sal4];
dat = [dat1;dat2;dat3;dat4];

died_r = dat(:,1)./(dat(:,1)+dat(:,2));


sal2 = 0:35;
S_ref = 24;
ALPHA = 6;
KSAL = 1;
MOR_SAL = (1-tanh((sal2/S_ref).^ALPHA))*KSAL;


figure;
scatter(sal,died_r,'r+');
hold on;
plot(sal2,MOR_SAL);