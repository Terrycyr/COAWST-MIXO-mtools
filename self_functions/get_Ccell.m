function [m_Ccell,a_Ccell] = get_Ccell
%Summary of this function goes here
%   Detailed explanation goes here


allo_a = 0.41;
allo_b = 0.7646;

m_ESD = 22; % um
%mixotroph cell radius (m)
m_rad=m_ESD/2./1000000.;
m_Vcell = 4./3.*3.14159*(m_ESD/2.)^3;
m_Ccell = allo_a/1000.*m_Vcell^allo_b*1e-6;

a_ESD = 1; % um
%prey cell radius (m)
a_rad=a_ESD/2./1000000.;
a_Vcell = 4./3.*3.14159*(a_ESD/2.)^3;
a_Ccell = allo_a/1000.*a_Vcell^allo_b*1e-6;

end