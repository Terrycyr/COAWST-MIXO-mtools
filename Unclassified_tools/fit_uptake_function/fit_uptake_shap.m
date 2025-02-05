clear all; close all;
%Kq
NC = 0.1:0.01:0.3;
Kq = 0.1;
Hq = 4;
NCabs = 0.3;
RF1 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

NC = 0.1:0.01:0.3;
Kq = 0.2;
Hq = 4;
NCabs = 0.3;
RF2 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

NC = 0.1:0.01:0.3;
Kq = 0.4;
Hq = 4;
NCabs = 0.3;
RF3 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

figure;
plot(NC,RF1);
hold on;
plot(NC,RF2);
hold on;
plot(NC,RF3);
hold on;
legend('Kq=0.1','Kq=0.2','Kq=0.4');

%Hq
NC = 0.1:0.01:0.3;
Kq = 0.3;
Hq = 1;
NCabs = 0.3;
RF1 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

NC = 0.1:0.01:0.3;
Kq = 0.4;
Hq = 4;
NCabs = 0.3;
RF2 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

NC = 0.1:0.01:0.3;
Kq = 0.4;
Hq = 8;
NCabs = 0.4;
RF3 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

figure;
plot(NC,RF1);
hold on;
plot(NC,RF2);
hold on;
plot(NC,RF3);
hold on;
legend('Hq=1','Hq=4','Hq=8');

%New&Original
NC = 0.1:0.01:0.3;
Kq = 0.1;
Hq = 4;
NCabs = 0.3;
RF1 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

NC = 0.1:0.01:0.3;
Kq = 0.4;
Hq = 8;
NCabs = 0.3;
RF2 = (((1+Kq^Hq)*(1-NC/NCabs).^Hq))./((1-NC/NCabs).^Hq+Kq^Hq);

figure;
plot(NC,RF1);
hold on;
plot(NC,RF2);
hold on;
legend('Original','New');