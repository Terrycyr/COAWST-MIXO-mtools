clear all; close all;


Rnutri = 0:0.01:1;
KNUTR = 0.6;
MOR_NUTR1 = (1-Rnutri.^2)*KNUTR;

Rnutri = 0:0.01:1;
R_ref = 0.22;
ALPHA = 2;
KNUTR = 0.9;
MOR_NUTR2 = (1-tanh((Rnutri/R_ref).^ALPHA))*KNUTR;

Rnutri = 0:0.01:1;
R_ref = 0.25;
ALPHA = 1;
KNUTR = 0.9;
MOR_NUTR3 = (1-tanh((Rnutri/R_ref).^ALPHA))*KNUTR;


figure;
plot(Rnutri,MOR_NUTR1);
hold on;
plot(Rnutri,MOR_NUTR2);
hold on;
plot(Rnutri,MOR_NUTR3);
legend('Old','New','New2');