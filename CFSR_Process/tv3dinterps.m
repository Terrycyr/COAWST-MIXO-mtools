function [ out_dat ] = tv3dinterps( data, input_x, input_y, out_x, out_y,method)
%v3dinterps Summary of this function goes here
%   Detailed explanation goes here
    print_int = linspace(0,100,41);
    out_dat = zeros([size(out_x) size(data,3)]);
    for t = 1:size(data,3)
        if(sum(t*100/size(data,3)>print_int))
               disp(['Spatial Interpolation Progress:', sprintf('%6.2f',t*100/size(data,3)),'%'])
               [~,pos] = max((t*100/size(data,3))-print_int);
               print_int(pos) = 1000;
        end
        
        tmp2 = data(:,:,t);
        F1 = TriScatteredInterp(input_x(:), input_y(:), tmp2(:),method);
        tmp = F1(out_x,out_y); 
        F2 = TriScatteredInterp(input_x(~isnan(tmp2)), input_y(~isnan(tmp2)),tmp2(~isnan(tmp2)),'nearest');
        tmp(isnan(tmp)) =  F2(out_x(isnan(tmp)), out_y(isnan(tmp)));
        
        out_dat(:,:,t) = tmp;
    end
end

