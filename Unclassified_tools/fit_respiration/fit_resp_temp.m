clear all; close all;

t = 0:35;

resp1 = 0.02*1.15.^(t-20);
resp2 = 0.015*1.15.^(t-20);
resp3 = 0.05*1.08.^(t-20);

figure;
plot(t,resp1);
hold on;
plot(t,resp2);
hold on;
plot(t,resp3);
legend('Kb','Kb2','Syn'); 