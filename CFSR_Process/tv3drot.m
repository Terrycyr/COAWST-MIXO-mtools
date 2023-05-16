function [ out_u, out_v ] = tv3drot( data_u, data_v , angle )
%tv3drot Summary of this function goes here
%   Detailed explanation goes here
%   angle, the angle between eastward and the data x-axis, unit: rad 
    out_u = zeros(size(data_u));
    out_v = zeros(size(data_v));
    for t = 1:size(data_u,3)
        out_u(:,:,t) = cos(angle).*data_u(:,:,t)+sin(angle).*data_v(:,:,t);
        out_v(:,:,t) = -sin(angle).*data_u(:,:,t)+cos(angle).*data_v(:,:,t);
    end       
end

