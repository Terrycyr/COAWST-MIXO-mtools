function [el,t,s_lon,s_lat,s_name] = read_noaa_tide(fname)
%read_noaa_tide Summary of this function goes here
%   Detailed explanation goes here
    [dat,txt,raw]  =xlsread(fname);
    s_lat = dat(1,6)+dat(1,7)/60;
    s_lon = -1*(dat(2,6)+dat(2,7)/60);
    el = dat(:,4);
    for i=2:length(el)+1
        tmp = strsplit(txt{i,1},'/');
        year(i) = str2num(tmp{3});
        mon(i) = str2num(tmp{1});
        day(i) = str2num(tmp{2});
        t(i-1) = datenum(year(i),mon(i),day(i));
        t(i-1) = t(i-1)+dat(i-1,1);
    end
    s_name = raw{1,7};   
end

