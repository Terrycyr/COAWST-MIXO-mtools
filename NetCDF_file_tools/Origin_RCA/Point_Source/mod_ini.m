clear all;

nh4t = ncread("Y2000_eutr_ini.nc",'NH4T');

nh4t = nh4t*0.0;

ncwrite("Y2000_eutr_ini.nc",'NH4T',nh4t);

nh4t = ncread("Y2000_eutr_ini.nc",'NH4T');