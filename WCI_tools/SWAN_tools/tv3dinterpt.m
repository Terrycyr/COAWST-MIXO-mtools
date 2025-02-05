function [ out ] = tv3dinterpt( data , origin_date, out_date, method )
%v3dinterpt Summary of this function goes here
%   Detailed explanation goes here
   print_int = linspace(0,100,41);
   [r,c,t] = size(data);
   out = zeros(r,c,length(out_date));
   k = 0;
   for i = 1:r
       for j = 1:c
           k = k+1;
           if(sum(k*100/(r*c)>print_int))
               disp(['Temporal Interpolation Progress:', sprintf('%6.2f',k*100/(r*c)),'%']);
               [~,pos] = max((k*100/(r*c))-print_int);
               print_int(pos) = 1000;
           end
           origin_t = origin_date;
           origin_d = squeeze(data(i,j,:));
           tmp = interp1(origin_t, origin_d, out_date, method);
           tmp(isnan(tmp)) = interp1(origin_t(~isnan(origin_d)), origin_d(~isnan(origin_d))...
               , out_date(isnan(tmp)), 'nearest','extrap');
           out(i,j,:) = tmp;
       end
   end
end

