function [ Dkm ] = distance( Lat1,Lon1,Lat2,Lon2 )
%distance Summary of this function goes here
%   Detailed explanation goes here
    p = pi/180;
    c1 = cos((Lat2 - Lat1) * p);
    c2 = cos(Lat1 * p);
    c3 = cos(Lat2 * p);
    c4 = cos((Lon2 - Lon1) * p);
    a = 0.5 - c1./2 + c2 .* c3 .* (1 - c4) ./ 2;
    Dkm = 12742 * asin(sqrt(a)); % in km
end

