clear all;

t = 0:35;

%Lenes et al., 2012
GITMAX1 = 1.8*exp(0.0633*(t-27));

TOPT2 = 30;
K2C = 1.8;
K2BETA1 = 0.003;
K2BETA2 = 0.001;
for i=1:length(t)
    if(t(i)<=TOPT2)
        GITMAX2(i) = K2C*exp(-K2BETA1*(t(i)-TOPT2)^2.);
    else
        GITMAX2(i) = K2C*exp(-K2BETA2*(TOPT2-t(i))^2.);
    end
end

figure;
plot(t,GITMAX1,t,GITMAX2);
legend('Lenes','Fitted');