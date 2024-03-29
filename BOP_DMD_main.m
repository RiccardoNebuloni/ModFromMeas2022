clear;close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simple example --- how to call code
%
% Here we fit data generated from 3 
% spatial modes, each with time dynamics 
% which are exponential in time
%
% The examples show how to call the optdmd
% wrapper with various options
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate synthetic data

% set up modes in space

x0 = 0;
x1 = 1;
nx = 200;

% space

xspace = linspace(x0,x1,nx);

% modes

f1 = sin(xspace);
f2 = cos(xspace);
f3 = tanh(xspace);

% set up time dynamics

t0 = 0;
t1 = 1;
nt = 100;

ts = linspace(t0,t1,nt);

% eigenvalues

e_1 = 1;
e_2 = -2;
e_3 = 1i;

evals = [e_1;e_2;e_3];

% create clean dynamics

xclean = f1'*exp(e_1*ts) + f2'*exp(e_2*ts) + f3'*exp(e_3*ts);

% %corrupt results
% sam= 60;
% Atest2 = zeros(nx,length(ts));
% Arand1 = rand(nx,length(ts));
% Arand2 = rand(nx,length(ts));
% r1k = randperm(length(ts)*n,sam);
% 
% for j = 1:sam
%         Atest2(r1k(j))=1;
% end
% Anoise = Atest2.*(i*Arand1);
% Xnoise = xclean+Anoise;



%% try bootstrapping - need to decide with or without replacement
%set seed
rng(7);


%number of time points
n = length(ts);

%number you want to choose
p = 50;

%number of noise cycles
num_Noisecycles = 1;

%number of cycles for each noise cycle
num_cycles =  100;

r = 3;


%create lambda vec for DMD cycles
lambda_vec_DMD = zeros(r,num_Noisecycles);

%create lambda vec for optdmd cycles
lambda_vec_optDMD = zeros(r,num_Noisecycles);



%create lambda vec for optdmd cycles
lambda_vec_mean_ensembleDMD = zeros(r,num_Noisecycles);

for k = 1: num_Noisecycles

%create data for noise cycle

sigma = .05;
xdata = xclean + sigma*randn(size(xclean));

%create lambda vec for ensembleDMD cycle
lambda_vec_ensembleDMD = zeros(r,num_cycles);
b_vec_ensembleDMD = zeros(r,num_cycles);
w_vec_ensembleDMD = zeros(length(xspace),r,num_cycles);

%try DMD
% [phi_DMD, lam_DMD, b_DMD, sig_DMD]= DMD(xdata(:,1:end-1), xdata(:,2:end), 3);


[phi_DMD, lam_DMD, b_DMD]= DMD(xdata(:,1:end-1), xdata(:,2:end), r);

%try regular optdmd
[w_opt,e_opt,b_opt] = optdmd(xdata,ts,r,1, varpro_opts('ifprint',0));
x1_opt = w_opt*diag(b_opt)*exp(e_opt*ts);

figure(1)
hold on
grid on
waterfall(xspace,ts',abs(xdata.'))
hold off

figure(2)
hold on
grid on
waterfall(xspace,ts',abs(x1_opt.'))
hold off


for j = 1:num_cycles
        %try with ioptdmd with DMD modes/evals as IC
        %select indices
        unsorted_ind = randperm(n,p);
        %sort ind so in ascending order. NOTE: evals have variable delta t
        ind = sort(unsorted_ind);

        %create dataset for this cycle by taking aforementioned indices
        xdata_cycle = xdata(:,ind);
        %selected index times
        ts_ind = ts(ind);

        [w_cycle,e1_cycle,b_cycle] = optdmd(xdata_cycle,ts_ind,r,1,varpro_opts('ifprint',0),e_opt);
        lambda_vec_ensembleDMD(:,j) = e1_cycle;
        b_vec_ensembleDMD(:,j) = b_cycle;
        w_vec_ensembleDMD(:,:,j) = w_cycle;
end

%store values
lambda_vec_DMD(:,k) = diag(lam_DMD);
lambda_vec_optDMD(:,k) = e_opt;

lambda_average = mean(lambda_vec_ensembleDMD,2);
b_average = mean(b_vec_ensembleDMD,2);
w_average = mean(w_vec_ensembleDMD,3);

x1_BoP = w_average*diag(b_average)*exp(lambda_average*ts);

figure(3)
hold on
grid on
waterfall(xspace,ts',abs(x1_BoP.'))
hold off



% figure()
% set(groot, 'defaultLineMarkerSize',15)
% scatter(real(evals),imag(evals),'bo')
% hold on
% scatter(real(e1_cycle),imag(e1_cycle),'rd')
% legend({'true eigenvalues','opt dmd boost'},'Location','NorthWest')





%make dist for DMD
sortedLambda_DMD = sort(lambda_vec_DMD,1,'ComparisonMethod','real');

%make dist for optDMD
sortedLambda_optDMD = sort(lambda_vec_optDMD,1,'ComparisonMethod','real');



%make dist for DMD
sortedLambda_ensembleDMD = sort(lambda_vec_ensembleDMD,1,'ComparisonMethod','real');

%got rid of absolute value
mean_value_ensemble_1 = mean(sortedLambda_ensembleDMD(1,:));
mean_value_ensemble_2 = mean(sortedLambda_ensembleDMD(2,:));
mean_value_ensemble_3 = mean(sortedLambda_ensembleDMD(3,:));

figure(7)
histfit(abs(sortedLambda_ensembleDMD(1,:)),40);
hold on;
xline(abs(sortedLambda_optDMD(1,k)),'k:','linewidth', 4);
hold on ; 
xline(abs(mean_value_ensemble_1),'b--', 'linewidth',4);
hold on; 
xline(2,'r-', 'linewidth',4);
% title('EnsembleDMD eigenvalue 1')
ax= gca;
ax.FontSize = 20;

figure(8)
histfit(abs(sortedLambda_ensembleDMD(2,:)),40);
hold on;
xline(abs(sortedLambda_optDMD(2,k)),'k:','linewidth', 4);
hold on ; 
xline(abs(mean_value_ensemble_2),'b--', 'linewidth',4);
hold on; 
xline(1,'r-', 'linewidth',4);
% title('EnsembleDMD eigenvalue 2')
ax= gca;
ax.FontSize = 20;


figure(9)
histfit(real(sortedLambda_ensembleDMD(3,:)),40);
hold on;
opt = xline(abs(sortedLambda_optDMD(3,k)),'k:','linewidth', 4);
hold on ; 
bop = xline(abs(mean_value_ensemble_3),'b--', 'linewidth',4);
hold on; 
tru = xline(1,'r-', 'linewidth',4);
% title('EnsembleDMD eigenvalue 3')
ax= gca;
ax.FontSize = 20;
% legend([bop,opt, tru],'BOP-DMD','OptDMD','True Value')

lambda_vec_mean_ensembleDMD(:, k)  = [mean_value_ensemble_1; mean_value_ensemble_2; mean_value_ensemble_3];

if mod(k,10) ==0
    disp(['This is the noise cycle ' num2str(k)])
end 

end

% save('evecs_high_03_noabs.mat','lambda_vec_mean_ensembleDMD','sortedLambda_optDMD')

% exit;