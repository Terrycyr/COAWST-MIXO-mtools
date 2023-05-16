function [ out ] = tv2dinterpt( data , origin_date, out_date, method, extrap )
%v3dinterpt Summary of this function goes here
%   Detailed explanation goes here
   print_int = linspace(0,100,41);
   [r,~] = size(data);
   out = zeros(r,length(out_date));
   k = 0;
   for i = 1:r
       k = k+1;
       if(sum(k*100/(r)>print_int))
           disp(['Temporal Interpolation Progress:', sprintf('%6.2f',k*100/(r)),'%']);
           [~,pos] = max((k*100/(r))-print_int);
           print_int(pos) = 1000;
       end
       origin_t = origin_date;
       origin_d = squeeze(data(i,:));
       tmp = interp1(origin_t, origin_d, out_date, method);
       if(extrap)
           tmp(isnan(tmp)) = interp1(origin_t(~isnan(origin_d)), origin_d(~isnan(origin_d))...
               , out_date(isnan(tmp)), 'nearest','extrap');
       end
       out(i,:) = tmp;
   end
end

