clear all; close all;

salt = 0:0.1:35;

cff1 = 5.922E-05;
cff2 = -0.007961;
cff3 = 0.3451;
cff4  = -3.866;

RSAL1 = cff1*salt.^3+ cff2*salt.^2+ cff3*salt+ cff4;

for i=1:length(salt)
    if(salt(i)<=35.)
        RSAL2(i) = exp(-0.008*(salt(i)-35.).^2.);
    end
    if(salt(i)>35.)
        RSAL2(i) = exp(-0.001*(35-salt(i)).^2.);
    end
end

RSAL1 = max(RSAL1,0.);
RSAL1 = min(RSAL1,1.);

RSAL2 = max(RSAL2,0.);
RSAL2 = min(RSAL2,1.);

figure;
plot(salt,RSAL1);
hold on;
plot(salt,RSAL2);