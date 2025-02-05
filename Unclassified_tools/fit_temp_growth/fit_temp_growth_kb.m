clear all; close all

dat = dlmread('tvsgrowth.txt',',');

t = 10:35;
% 
%Lenes et al., 2012
GITMAX1 = 0.8*exp(0.0633*(t-27));

TOPT2 = 27;
K2C = 0.8;
K2BETA1 = 0.01;
K2BETA2 = 0.03;
for i=1:length(t)
    if(t(i)<=TOPT2)
        GITMAX2(i) = K2C*exp(-K2BETA1*(t(i)-TOPT2)^2.);
    else
        GITMAX2(i) = K2C*exp(-K2BETA2*(TOPT2-t(i))^2.);
    end
end

figure;
scatter(dat(:,1),K2C*dat(:,2),'k','filled');
hold on;
plot(t,GITMAX1,t,GITMAX2);
set(gcf,'color','w');
legend('Tilney et al., 2019','Lenes et al., 2012','Model');