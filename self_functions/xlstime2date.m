function [date] = xlstime2date(time)
%xlstime2date Summary of this function goes here
%   Detailed explanation goes here
   date = time+datenum(1900,1,0)-1;
end