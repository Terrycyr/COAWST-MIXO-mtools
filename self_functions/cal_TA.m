function TA = cal_TA(pH,DIC,sal,tin,tout,pin,pout,si,sol_p,9)
%get_station_ij Summary of this function goes here
%   Detailed explanation goes here

addpath(path,'C:\Users\cheny\Desktop\matlab_tools\CO2SYS');

%      indicating that PAR1 (or PAR2) is of type:
%  1 = Total Alkalinity
%  2 = DIC
%  3 = pH
%  4 = pCO2
%  5 = fCO2

%      indicating the K1 K2 dissociation constants that are to be used:
%   1 = Roy, 1993											T:    0-45  S:  5-45. Total scale. Artificial seawater.
%   2 = Goyet & Poisson										T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
%   3 = HANSSON              refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   4 = MEHRBACH             refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO	T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   6 = GEOSECS (i.e., original Mehrbach)					T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   7 = Peng	(i.e., originam Mehrbach but without XXX)	T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)	T:    0-50  S:     0. 
%   9 = Cai and Wang, 1998									T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
%  10 = Lueker et al, 2000									T:    2-35  S: 19-43. Total scale. Real seawater.
%  11 = Mojica Prieto and Millero, 2002.					T:    0-45  S:  5-42. Seaw. scale. Real seawater
%  12 = Millero et al, 2002									T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
%  13 = Millero et al, 2006									T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  14 = Millero et al, 2010									T:    0-50  S:  1-50. Seaw. scale. Real seawater.

par1type = ones(size(DIC))*3;
par1    =   pH;    % value of the first parameter
par2type = ones(size(DIC))*2;
par2    =   DIC;    % value of the second parameter
sal     =  sal;    % Salinity of the sample
tempin  =  tin;    % Temperature at input conditions
tempout = tout;    % Temperature at output conditions
presin  =  pin;    % Pressure    at input conditions
presout = pout;    % Pressure    at output conditions
sil     =   si;    % Concentration of silicate  in the sample (in umol/kg)
po4     =sol_p;    % Concentration of phosphate in the sample (in umol/kg)
pHscale = ones(size(DIC))*1;    % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c   = ones(size(DIC))*9;    % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c   = ones(size(DIC))*1;    % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")

% Do the calculation. See CO2SYS's help for syntax and output format
A=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,presin,presout,sil,po4,pHscale,k1k2c,kso4c);

TA = A(:,1); %umol/kg

end