clear all; close all


m_NC = 0.05:0.001:0.3;
m_NCabs = 0.3;
m_NCmin = 0.1;
Minval = 1e-7;

TERM = ((m_NCabs-m_NC)/(m_NCabs-m_NCmin+Minval))...
    ./((m_NCabs-m_NC)/(m_NCabs-m_NCmin+Minval)+0.01);

figure;
plot(m_NC,TERM);

ylim([-4 4]);