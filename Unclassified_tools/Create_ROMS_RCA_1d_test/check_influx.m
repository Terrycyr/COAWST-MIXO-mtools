clear all;

flow = ncread('./bio1dtest_river_bio.nc','river_transport');

shape = ncread('bio1dtest_river_bio.nc','river_Vshape');

influx = flow(end,:);

influx2 = influx*3600;

for i  = 1:length(influx2)
    influx_ac(i) = sum(influx2(1:i));
end

figure;
plot(influx_ac);