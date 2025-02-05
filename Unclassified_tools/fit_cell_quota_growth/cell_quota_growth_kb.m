clear all;close all;

%N:C, Corcoran et al., 2014
dat1 = [0.1 1/15;
       0.21 1/13
       0.23 1/7.2];

%Hardison et al., 2012
dat2 = [0.436 0.13;
        0.094 0.079;
        0.319 0.14;
        0.095 0.061];

%Sophia thesis, show a range of N:C (mass ratio) of 0.0778-0.134

dat = [dat1;dat2];
dat(:,1) = dat(:,1)/max(dat(:,1));
x = dat(:,2)*14/12;
y = dat(:,1);

m_KQN = 10;
m_NC = 0.09:0.001:0.17;
m_NCmax = 0.17;
m_NCmin = 0.09;
Minval = 1e-7;
m_NCu=(1+m_KQN)*(m_NC-m_NCmin)...
    ./((m_NC-m_NCmin)+m_KQN*(m_NCmax-m_NCmin)+Minval);

m_KQN2 = 10;
m_NC2 = 0.15:0.001:0.25;
m_NCmax2 = 0.25;
m_NCmin2 = 0.15;
m_NCu2=(1+m_KQN2)*(m_NC2-m_NCmin2)...
    ./((m_NC2-m_NCmin2)+m_KQN2*(m_NCmax2-m_NCmin2)+Minval);

figure;
plot(m_NC,m_NCu);
hold on;
%plot(m_NC2,m_NCu2);
%hold on;
scatter(x,y,'g','filled');

%%
clear all;
%P:C, Corcoran et al., 2014
dat = [0.1 290;
       0.21 210;
       0.23 140];
dat(:,1) = dat(:,1)/max(dat(:,1));
x = (1./dat(:,2))*31/12;
y = dat(:,1);

m_KQP = 0.1;
m_PC = 0.01:0.001:0.02;
m_PCmax = 0.02;
m_PCmin = 0.01;
Minval = 1e-7;


m_PCu=(1+m_KQP)*(m_PC-m_PCmin)...
    ./((m_PC-m_PCmin)+m_KQP*(m_PCmax-m_PCmin)+Minval);


figure;
plot(m_PC,m_PCu);
hold on;
scatter(x,y,'g','filled');
