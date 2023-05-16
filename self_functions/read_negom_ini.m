function [dat_final] = read_negom_ini( fn, lon, lat, sc, h, sheet)
%   read_negom_ini Summary of this function goes here
%   Detailed explanation goes here
    NEGOM = xlsread(fn,sheet);
    dep = NEGOM(:,2);
    dat_m = NEGOM(:,1);

    figure;
    plot(dep,dat_m);
    
    [r,c] = size(lon);
    [~,~,N] = size(sc);
    for i=1:r
        for j=1:c
            for k=1:N
                dat_final(i,j,k) = interp1(dep, dat_m,h(i,j).*sc(i,j,k),'linear')/1.025; %umol/l --> umol/kg
            end
        end
    end
end