clear all

t1 = 25;
t2 = 32;

g1 = 0.10*1.08^(t1-20);
g2 = 0.30*1.08^(t2-20);

x = [t1 t2];
y = [g1 g2];

t = 0:35;

graz1 = 0.1*1.08.^(t-20);
graz2 = 0.3*1.08.^(t-20);
%graz3 = 0.07*1.22.^(t-20);
graz3 = 0.08*1.15.^(t-20);

figure;
plot(t,graz1);
hold on;
plot(t,graz2);
hold on;
plot(t,graz3);
hold on;
scatter(x,y,'black','filled');
legend('2021','2022','New','2021&2022');