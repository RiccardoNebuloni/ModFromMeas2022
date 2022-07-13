%basata sul video
clear;close all;clc;
%% optdmd
addpath('./src');
addpath('./utils/');
load("input_data.mat"); %year, snowshoe hare pelts, lynx pelts
inputdata = inputdata';

space = [1 2]; %lepri e linci

% t = inputdata(1,1:end);

t_plot=inputdata(1,1:end);
dt = t_plot(2)-t_plot(1);
t = (t_plot-t_plot(1))/dt;
% t = t_plot;
pelts = inputdata(2:end,:);
[size_1,size_2] = size(pelts);

%figure(1), waterfall(space,t',pelts');
%figure(2), plot(t,pelts(1,:)');



X = pelts; param = 1;
X1 = X(:,1:end-param);
X2 = X(:,param+1:end);


r = 2;

% 1 --- fit to unprojected data

imode = 1;
[w,e1,b] = optdmd(X,t,r,imode);

x1_opt = w*diag(b)*exp(e1*t);
figure (1)
hold on
set(gca,'Fontsize',20)
grid on
plot(t_plot,1000*pelts(1,:)','b-o','LineWidth',2)
plot(t_plot',1000*abs(x1_opt(1,:)'),'c--','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Snowshore Hare')
legend('Real Data','Opt. DMD Approx.')
hold off

figure (2)
hold on
grid on
set(gca,'Fontsize',20)
plot(t_plot,1000*pelts(2,:)','r-x','LineWidth',2)
plot(t_plot',1000*abs(x1_opt(2,:)'),"Color",[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Lynx')
legend('Real Data','Opt. DMD Approx.')
hold off
%% bagging

num_cycles =  100;
%number of time points
n = length(t);

%number you want to choose
px = 2;
rx = 2;
lambda_vec_ensembleDMD = zeros(r,num_cycles);
b_vec_ensembleDMD = zeros(r,num_cycles);
w_vec_ensembleDMD = zeros(length(space),r,num_cycles);

for j = 1:num_cycles 
        %try with ioptdmd with DMD modes/evals as IC
        %select indices
        unsorted_ind = randperm(n,px);
        %sort ind so in ascending order. NOTE: evals have variable delta t
        ind = sort(unsorted_ind);

        %create dataset for this cycle by taking aforementioned indices
        xdata_cycle = X(:,ind);
        %selected index times
        t_ind = t(ind);

        [w_cycle,e1_cycle,b_cycle] = optdmd(xdata_cycle,t_ind,rx,1,varpro_opts('ifprint',0),e1);
%         [w_cycle,e1_cycle,b_cycle] = optdmd(xdata_cycle,t_ind,rx,imode);
        lambda_vec_ensembleDMD(:,j) = e1_cycle;
        b_vec_ensembleDMD(:,j) = b_cycle;
        w_vec_ensembleDMD(:,:,j) = w_cycle;
end

sortedLambda_ensembleDMD = sort(lambda_vec_ensembleDMD,1,'ComparisonMethod','real');
lambda_average = mean(lambda_vec_ensembleDMD,2);
b_average = mean(b_vec_ensembleDMD,2);
w_average = mean(w_vec_ensembleDMD,3);

x1_BoP = w_average*diag(b_average)*exp(lambda_average*t);
x1_BoP = real(x1_BoP);
figure (3)
hold on
grid on
plot(t_plot,pelts(1,:)')
plot(t_plot',abs(x1_BoP(1,:)'))
hold off

figure (4)
hold on
grid on
plot(t_plot,pelts(2,:)')
plot(t_plot',abs(x1_BoP(2,:)'))
hold off

figure (3)
hold on
set(gca,'Fontsize',20)
grid on
plot(t_plot,1000*pelts(1,:)','b-o','LineWidth',2)
plot(t_plot',1000*abs(x1_BoP(1,:)'),'c--','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Snowshore Hare')
legend('Real Data','Bagging Opt. DMD')
hold off

figure (4)
hold on
grid on
set(gca,'Fontsize',20)
plot(t_plot,1000*pelts(2,:)','r-x','LineWidth',2)
plot(t_plot',1000*abs(x1_BoP(2,:)'),"Color",[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Lynx')
legend('Real Data','Bagging Opt. DMD')
hold off

%% Time-delay model:
p = 6; contatore=0; %p=6
for h_idx = 1:p
    H(h_idx+contatore:h_idx*size(X,1),:) = X(:,1+contatore:end-p+contatore);
    contatore = contatore+1;
end

rH = 12; rH=min(rH,size(H,1));

t_plotH=t_plot(1:size(H,2));
tH = (t_plotH-t_plotH(1))/dt;
[UH,SH,VH] = svd(H,"econ");

figure
hold on
grid on
plot(sum(SH)./sum(sum(SH))*100,'k-o')
hold off

[wH,e1H,bH] = optdmd(H,tH,rH,imode);

% x1_optH = wH*diag(bH)*exp(e1H*tH);
x1_optH = wH*diag(bH)*exp(e1H*t); %tutto l'intervallo

figure (5)
hold on
set(gca,'Fontsize',20)
grid on
% plot(t_plotH,1000*pelts(1,1:size(H,2))','b-o','LineWidth',2)
plot(t_plot,1000*pelts(1,:)','b-o','LineWidth',2)
% plot(t_plotH',1000*abs(x1_optH(1,:)'),'c--','LineWidth',2)
plot(t_plot,1000*abs(x1_optH(1,:)'),'c--','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Snowshore Hare')
legend('Real Data','Time-Delay DMD')
hold off

figure (6)
hold on
set(gca,'Fontsize',20)
grid on
plot(t_plot,1000*pelts(2,:)','r-x','LineWidth',2)
plot(t_plot,1000*abs(x1_optH(2,:)'),"Color",[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Lynx')
legend('Real Data','Time-Delay DMD')
hold off
 
figure (7)
hold on
set(gca,'Fontsize',20)
grid on
plot(t_plot,1000*pelts(1,:)','b-o','LineWidth',2)
plot(t_plot,1000*abs(x1_optH(1,:)'),'c--','LineWidth',2)
plot(t_plot,1000*pelts(2,:)','r-x','LineWidth',2)
plot(t_plot,1000*abs(x1_optH(2,:)'),"Color",[0.4660 0.6740 0.1880],'LineStyle','-.','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Real Data vs Time-Delay DMD - d=6, r=12')
legend('Snowshoe Hare: Real Data','Snowshoe Hare: Time-Delay DMD','Lynx: Real Data','Lynx: Time-Delay DMD')
hold off
%% bagging time delay DMD

num_cyclesH =  600; %300
%number of time points
nH = length(tH);

%number you want to choose
pxH = 24;
lambda_vec_ensembleDMDH = zeros(rH,num_cyclesH);
b_vec_ensembleDMDH = zeros(rH,num_cyclesH);
w_vec_ensembleDMDH = zeros(size(H,1),rH,num_cyclesH);

for jH = 1:num_cyclesH 
        %try with ioptdmd with DMD modes/evals as IC
        %select indices
        unsorted_indH = randperm(nH,pxH);
        %sort ind so in ascending order. NOTE: evals have variable delta t
        indH = sort(unsorted_indH);

        %create dataset for this cycle by taking aforementioned indices
        xdata_cycleH = H(:,indH);
        %selected index times
        t_indH = tH(indH);

        [w_cycleH,e1_cycleH,b_cycleH] = optdmd(xdata_cycleH,t_indH,rH,1,varpro_opts('ifprint',0),e1H);
        lambda_vec_ensembleDMDH(:,jH) = e1_cycleH;
        b_vec_ensembleDMDH(:,jH) = b_cycleH;
        w_vec_ensembleDMDH(:,:,jH) = w_cycleH;
end

sortedLambda_ensembleDMDH = sort(lambda_vec_ensembleDMDH,1,'ComparisonMethod','real');
lambda_averageH = mean(lambda_vec_ensembleDMDH,2);
b_averageH = mean(b_vec_ensembleDMDH,2);
w_averageH = mean(w_vec_ensembleDMDH,3);

% x1_BoPH = w_averageH*diag(b_averageH)*exp(lambda_averageH*tH);
x1_BoPH = w_averageH*diag(b_averageH)*exp(lambda_averageH*t);


figure (8)
hold on
grid on
% plot(t_plotH,pelts(1,1:size(H,2))')
% plot(t_plotH',abs(x1_BoPH(1,:)'))
plot(t_plot,pelts(1,:)')
plot(t_plot',abs(x1_BoPH(1,:)'))
legend('real Data','Fitting')
hold off
% 
figure (9)
hold on
grid on
plot(t_plot,pelts(2,:)')
% plot(t_plot',abs(x1_BoPH(2,:)'))
plot(t_plot',XXX)
legend('real Data','Fitting')
hold off

figure(10)
hold on
set(gca,'Fontsize',20)
grid on
plot(t_plot,1000*pelts(1,:)','b-o','LineWidth',2)
plot(t_plot',1000*abs(x1_BoPH(1,:)'),'c--','LineWidth',2)
plot(t_plot,1000*pelts(2,:)','r-x','LineWidth',2)
plot(t_plot',1000*XXX,"Color",[0.4660 0.6740 0.1880],'LineStyle','-.','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Real Data vs Time-Delay BOP-DMD')
legend('Snowshoe Hare: Real Data','Snowshoe Hare: Time-Delay BOP-DMD','Lynx: Real Data','Lynx: Time-Delay BOP-DMD')
hold off
%% Lotka-Volterra Model fitting
%done
%x1 = (a1-a2y)x b=0.6156 p=0.0301
%y1 = (b2x-b1)y r=0.5706 d=0.0114
% a1 = 0.618; a2 = 0.0289; b2 = 0.0118; b1 = 0.5762;
a1 = 0.85; a2 = 0.0289; b2 = 0.0118; b1 = 0.41;

shoe = X(1,:)'; %x
linx = X(2,:)'; %y

%derivate
der_shoe = (repmat(a1,size(shoe,1),1) - repmat(a2,size(shoe,1),1).*linx).*shoe;
der_lynx = (repmat(b2,size(shoe,1),1).*shoe - repmat(b1,size(shoe,1),1)).*linx;
dx = [der_shoe der_lynx];


%forse va fatto sui dati veri.
dt=1;
for tt =1:size(X,2)-1
der_shoe_data(tt) = (X(1,tt+1)-X(1,tt))/dt; 
der_lynx_data(tt) = (X(2,tt+1)-X(2,tt))/dt; 
end

dx_data = [der_shoe_data.' der_lynx_data.'];
%% SINDY
polyorder = 2;  % search space up to fifth order polynomials
usesine = 0;    % no trig functions

n = 2; %order of the system

%u = reshape(X(:,2:end-1).',(size_2-2)*size_1,1);

Theta = poolData(X(:,1:end-1)',n,polyorder,usesine);
%Theta(:,3) = [];

%regression
%xi_inv=Theta\dx_data;
xi1_inv=Theta\der_shoe_data.';
xi2_inv=Theta\der_lynx_data.';

xi1_pinv=pinv(Theta)*der_shoe_data.';
xi2_pinv=pinv(Theta)*der_lynx_data.';


xi1_lasso=lasso(Theta,der_shoe_data,'Lambda',0.1);
xi2_lasso=lasso(Theta,der_lynx_data,'Lambda',0.1);

% xi1_lasso=lasso(Theta,der_shoe(1:end-1),'Lambda',0.1);
% xi2_lasso=lasso(Theta,der_lynx(1:end-1),'Lambda',0.1);

% Xi = [xi1_inv xi2_inv];
Xi = [xi1_lasso xi2_lasso];


% compute Sparse regression: sequential least squares
 lambda = 0.002;      % lambda is our sparsification knob.
 %Xi = sparsifyDynamics(Theta,dx_data,lambda,n);
%poolDataLIST({'Shoe','Lynx'},Xi,n,polyorder,usesine);

%ricostruzione approxximata
tspan = [0 30]; options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
t(1) = 0;
t=t'; x_piccolo = [shoe linx];
%% Bagging SINDY
num_cyclesS =  300;
%number of time points
nS = length(tB);
%number you want to choose
% pxS = round(0.85*size(Theta,2)); 
pxS = 7;
ensT = 0.65;
mOutBS = zeros(pxS,n,num_cyclesS);
libOutBS = zeros(pxS,num_cyclesS);

for jS = 1:num_cyclesS
     rs = RandStream('mlfg6331_64','Seed',jS);
     libOutBS(:,jS) = datasample(rs,1:size(Theta,2),pxS,'Replace',false)';
    mOutBS(:,:,jS) = sparsifyDynamics(Theta(:,libOutBS(:,jS)),dx_data,lambda,n);
end

inclProbBS = zeros(size(Theta,2),n);
for iii = 1:num_cyclesS
    for jjj = 1:n
        for kkk = 1:pxS
            if mOutBS(kkk,jjj,iii) ~= 0
                inclProbBS(libOutBS(kkk,iii),jjj) = inclProbBS(libOutBS(kkk,iii),jjj) + 1;
            end
        end
    end
end
inclProbBS = inclProbBS/num_cyclesS*size(Theta,2)/pxS;

XiD = zeros(size(Theta,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;
%     XiBias = sparsifyDynamics(Theta(:,libEntry),dx_data(:,iii),lambda,1);
    XiBias=lasso(Theta(:,libEntry),dx_data(:,iii),'Lambda',0.1);
    XiD(libEntry,iii) = XiBias;
end