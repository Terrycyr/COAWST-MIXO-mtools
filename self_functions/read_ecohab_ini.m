function [dat_final] = read_ecohab_ini( fn, lon, lat, sc, h, sheet)
%   read_negom_ini Summary of this function goes here
%   Detailed explanation goes here
    ECOHAB = xlsread(fn,sheet);
    dep = ECOHAB(:,1);
    dat_m = ECOHAB(:,2);

    figure;
    plot(dep,dat_m);
    
    [r,c] = size(lon);
    [~,~,N] = size(sc);
    for i=1:r
        for j=1:c
            for k=1:N
                dat_final(i,j,k) = interp1(dep, dat_m,h(i,j),'linear')/1.025; %umol/l --> umol/kg
            end
        end
    end
end