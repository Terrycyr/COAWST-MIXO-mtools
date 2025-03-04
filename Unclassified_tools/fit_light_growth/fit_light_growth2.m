clear all;close all;
%Diatom
x1 =[7.583947125
19.76643235
19.81684264
18.52297869
21.22833058
34.93992775
33.66286722
62.06066037
62.06066037
58.15386339
230.1117428
232.8002913
402.7417593
402.8005713
701.0866217
699.8347663
];

y1 = [0.057061641
0.234533845
0.272509592
0.297798751
0.335830509
0.665233708
0.703181449
1.096185622
1.096185622
1.153065225
1.694667712
1.720040888
1.742613493
1.786918531
1.495743132
1.552678746
]; %S.tropicum, Hideki Kaeriyama, 2011, South Japan

x2 = [15 
    38 
    130 
    650 
    1500];
y2 = [0.1 
    0.32 
    1. 
    1.99 
    1.68];

%
y1 = y1/max(y1);
y2 = y2/max(y2);
%y3 = y3/max(y3);

x = [x1;x2];
y = [y1;y2];

L=0:1700;
Ls = 190;
mu1 = L./Ls.*exp(1-(L./Ls));

alpha2 = 0.02;
K2C = 2.3; 
mu2 = tanh(alpha2*L/K2C);

%(1-exp(-a.*x./c)).*exp(-b.*x./c)
alpha3 = 3.5;
beta3 = 0.09;
Ls = 275;
mu3 = (1-exp(-alpha3.*L./Ls)).*exp(-beta3.*L./Ls);
mu_max =(1-exp(-alpha3*(-Ls/alpha3*log(beta3/(beta3+alpha3)))/Ls))*exp(-beta3*(-Ls/alpha3*log(beta3/(beta3+alpha3)))/Ls);
mu3 = mu3./mu_max;

Linhib = 500;
alpha4 = 0.012;
beta4 = 0.0003;
for i=1:length(L)
    L_inhib = L(i)-Linhib;
    mu4(i) = tanh(alpha4*L(i));
    if(L_inhib>0)
        mu4(i) = mu4(i)*exp(-beta4*L_inhib);
    end
end



figure('Units','pixels','Position',[100 100 1500 400]);
t = tiledlayout(1,2);
nexttile;
scatter(x1,y1,'filled','r');
hold on;
scatter(x2,y2,'filled','g');
hold on;
%scatter(x3,y3,'filled','b','linewidth',2);
plot(L,mu1,L,mu2,L,mu3,L,mu4);

nexttile;
scatter(x1,y1,'filled','r');
hold on;
scatter(x2,y2,'filled','g');
hold on;
%scatter(x3,y3,'filled','b','linewidth',2);
plot(L,mu1,L,mu2,L,mu3,L,mu4);
xlim([0 200]);
legend('Data1','Data2','Walsh','Hyperbolic Tangent Function','CoSine','Hyperbolic Tangent Function+Light inhibition','location','southeast');
title(t,'Diatoms');

diatom_x = L;
diatom_y = mu4;

Linhib_all = Linhib;
alpha_all = alpha4;
beta_all = beta4;

% Synechococcus sp.
x1 = [6.198547215496369
5.230024213075062
21.210653753026634
25.084745762711876
21.45278450363196
55.35108958837773
55.108958837772406
46.150121065375316
70.12106537530266
99.66101694915255
120
];

y1 = [ 0.2108843537414966
0.4659863945578231
0.44557823129251695
0.7006802721088435
0.6938775510204082
0.8843537414965986
1
1.1598639455782314
1.207482993197279
1.0374149659863945
1.1462585034013606
]; %Lisa Campbell, 1986

x2 = [3.125
0.78125
0.546875
1.953125
6.875
7.8125
30.078125
31.484375
27.734375
21.875
47.1875
50.703125
55.859375
82.109375
53.28125
];

y2=[0.047311828
0.078853047
0.093189964
0.118996416
0.140501792
0.207885305
0.362724014
0.411469534
0.461648746
0.470250896
0.450179211
0.397132616
0.488888889
0.533333333
0.531899642
];    %Timmermans, 2005

x3 = [15 25 50 86 207 445 550 650]';

y3 = [0.52 0.61 0.72 0.94 1.13 1.04 0.86 0.89]'; %Six,2004

y1 = y1/max(y1);
y2 = y2/max(y2);
y3 = y3/max(y3);

x = [x1;x2;x3];
y = [y1;y2;y3];

L=0:1700;
Ls = 125;
mu1 = L./Ls.*exp(1-(L./Ls));

alpha2 = 0.02;
K2C = 1.0; 
mu2 = tanh(alpha2*L/K2C);

%(1-exp(-a.*x./c)).*exp(-b.*x./c)
alpha3 = 2.996;
beta3 = 0.02637;
Ls = 70;
mu3 = (1-exp(-alpha3.*L./Ls)).*exp(-beta3.*L./Ls);
mu_max =(1-exp(-alpha3*(-Ls/alpha3*log(beta3/(beta3+alpha3)))/Ls))*exp(-beta3*(-Ls/alpha3*log(beta3/(beta3+alpha3)))/Ls);
mu3 = mu3./mu_max;

Linhib = 200;
alpha4 = 0.02;
beta4 = 0.00045;
for i=1:length(L)
    L_inhib = L(i)-Linhib;
    mu4(i) = tanh(alpha4*L(i));
    if(L_inhib>0)
        mu4(i) = mu4(i)*exp(-beta4*L_inhib);
    end
end

figure('Units','pixels','Position',[100 100 1500 400]);
t = tiledlayout(1,2);
nexttile;
scatter(x1,y1,'filled','r');
hold on;
scatter(x2,y2,'filled','g');
hold on;
scatter(x3,y3,'filled','b','linewidth',2);
plot(L,mu1,L,mu2,L,mu3,L,mu4);

nexttile;
scatter(x1,y1,'filled','r');
hold on;
scatter(x2,y2,'filled','g');
hold on;
scatter(x3,y3,'filled','b','linewidth',2);
plot(L,mu1,L,mu2,L,mu3,L,mu4);
legend('Data1','Data1','Data3','Walsh','Hyperbolic Tangent Function','CoSine','Hyperbolic Tangent Function+Light inhibition','location','southeast');
title(t,'\itSyn.');
xlim([0 200]);

syn_x = L;
syn_y = mu4;


Linhib_all = [Linhib_all,Linhib];
alpha_all = [alpha_all,alpha4];
beta_all = [beta_all,beta4];

%K. brevis
x1 = [1339 1218 534 289 125 78 43]';
y1 = [0.49 0.52 0.52 0.53 0.36 0.24 0.14]';

x2 = [1339 1218 534 289 125 78 43 17 10]';
y2 = [0.5 0.45 0.42 0.48 0.37 0.33 0.22 0 0]';

x3 = [1339 1218 534 289 125 78 43 17 10]';
y3 = [0.47 0.48 0.49 0.47 0.36 0.3 0.11 0 0]';

raw_data4 = [18.70563674321503, -0.12334801762114533
18.70563674321503, -0.20264317180616742
18.70563674321503, -0.2819383259911894
18.70563674321503, -0.32422907488986785
30.396659707724424, -0.04933920704845807
30.396659707724424, -0.10220264317180605
30.396659707724424, -0.15506607929515415
51.44050104384134, 0.07753303964757718
51.44050104384134, 0.2044052863436125
51.44050104384134, 0.2784140969162997
66.63883089770356, 0.19383259911894285
66.93110647181628, 0.34185022026431733
66.63883089770356, -0.022907488986784075
122.75574112734863, 0.2044052863436125
122.75574112734863, 0.2678414096916301
122.4634655532359, 0.325991189427313
51.44050104384134, -0.20792951541850213];

x4 = raw_data4(:,1);
y4 = (raw_data4(:,2)-min(raw_data4(:,2)))/(max(raw_data4(:,2))-min(raw_data4(:,2)));

y1 = y1/max(y1);
y2 = y2/max(y2);
y3 = y3/max(y3);

x = [x1;x2;x3];
y = [y1;y2;y3];

L=0:1700;
Ls = 60;
mu1 = L./Ls.*exp(1-(L./Ls));

alpha2 = 0.02;
K2C = 0.8; 
mu2 = tanh(alpha2*L/K2C);

%(1-exp(-a.*x./c)).*exp(-b.*x./c)
alpha3 = 2.2;
beta3 = 0.002;
Ls = 60;
mu3 = (1-exp(-alpha3.*L./Ls)).*exp(-beta3.*L./Ls);
mu_max =(1-exp(-alpha3*(-Ls/alpha3*log(beta3/(beta3+alpha3)))/Ls))*exp(-beta3*(-Ls/alpha3*log(beta3/(beta3+alpha3)))/Ls);
mu3 = mu3./mu_max;


Linhib = 500;
alpha4 = 0.023;
beta4 = 0.00005;
for i=1:length(L)
    L_inhib = L(i)-Linhib;
    mu4(i) = tanh(alpha4*L(i));
    if(L_inhib>0)
        mu4(i) = mu4(i)*exp(-beta4*L_inhib);
    end
end

figure('Units','pixels','Position',[100 100 1500 400]);
t = tiledlayout(1,2);
nexttile
scatter(x1,y1,'filled','r');
hold on;
scatter(x2,y2,'filled','g');
hold on;
scatter(x3,y3,'filled','b','linewidth',2);
hold on;
scatter(x4,y4,'filled','y','linewidth',2);
plot(L,mu1,L,mu2,L,mu3,L,mu4);

nexttile
scatter(x1,y1,'filled','r');
hold on;
scatter(x2,y2,'filled','g');
hold on;
scatter(x3,y3,'filled','b','linewidth',2);
hold on;
scatter(x4,y4,'filled','y','linewidth',2);
plot(L,mu1,L,mu2,L,mu3,L,mu4);
legend('Data1-Tilney','Data2-Tilney','Data3-Tilney','Data4-Magana','Walsh','Hyperbolic Tangent Function','CoSine','Hyperbolic Tangent Function+Light inhibition','location','southeast');
title(t,'\itK.brevis');
xlim([0 200]);

kb_x = L;
kb_y = mu4;

Linhib_all = [Linhib_all,Linhib];
alpha_all = [alpha_all,alpha4];
beta_all = [beta_all,beta4];

%uflagellate
L=0:1700;
Ls = 275;
mu1 = L./Ls.*exp(1-(L./Ls));

K2C = 1.6;
Linhib = 300;
alpha4 = 0.0092;
beta4 = 0.0001;
for i=1:length(L)
    L_inhib = L(i)-Linhib;
    mu4(i) = tanh(alpha4*L(i));
    if(L_inhib>0)
        mu4(i) = mu4(i)*exp(-beta4*L_inhib);
    end
end
figure('Units','pixels','Position',[100 100 1000 400]);
plot(L,mu1,L,mu4);
legend('Walsh','Hyperbolic Tangent Function+Light inhibition','location','northeastoutside');
title('{\itu}flagellate');

u_x = L;
u_y = mu4;

Linhib_all = [Linhib_all,Linhib];
alpha_all = [alpha_all,alpha4];
beta_all = [beta_all,beta4];

figure;
plot(diatom_x,diatom_y,syn_x,syn_y,kb_x,kb_y,u_x,u_y,'LineWidth',2);
legend('Diatoms','Synechococcus sp.','K. brevis','{\itu}flagellate');
xlim([0 1700]);

disp(['    Diatom     ','Synechococcus     ','K. brevis     ','uflagellate'])
Linhib_all
alpha_all
beta_all
