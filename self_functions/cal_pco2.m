function pco2 = cal_pco2(TA,DIC,sal,tin,tout,pin,pout,si,sol_p)
%get_station_ij Summary of this function goes here
%   Detailed explanation goes here

addpath(path,'C:\Users\cheny\Desktop\matlab_tools\CO2SYS');

%      indicating that PAR1 (or PAR2) is of type:
%  1 = Total Alkalinity
%  2 = DIC
%  3 = pH
%  4 = pCO2
%  5 = fCO2

par1type =   ones(size(TA))*1;
par1    =   TA;    % value of the first parameter
par2type =   ones(size(TA))*2;
par2    = DIC;    % value of the second parameter
sal     =  sal;    % Salinity of the sample
tempin  =  tin;    % Temperature at input conditions
tempout = tout;    % Temperature at output conditions
presin  =  pin;    % Pressure    at input conditions
presout = pout;    % Pressure    at output conditions
sil     =   si;    % Concentration of silicate  in the sample (in umol/kg)
po4     =sol_p;    % Concentration of phosphate in the sample (in umol/kg)
pHscale =    ones(size(TA))*1;    % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c   =    ones(size(TA))*9;    % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("9" means "Cai and Wang, 1998")
kso4c   =    ones(size(TA))*1;    % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

pco2 = A(:,19); %atm
end