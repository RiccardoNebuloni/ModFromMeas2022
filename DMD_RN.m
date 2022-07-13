%basata sul video
clear;close all;clc;
%%
load("input_data.mat"); %year, snowshoe hare pelts, lynx pelts
inputdata = inputdata';

space = [1 2]; %lepri e linci

% t = inputdata(1,1:end);
t = 0:29;
pelts = inputdata(2:end,:);

%figure(1), waterfall(space,t',pelts');
%figure(2), plot(t,pelts(1,:)');


dt = t(2)-t(1);
X = pelts; param = 1;
X1 = X(:,1:end-param);
X2 = X(:,param+1:end);
u = X(:,1); %initial condition
[U,Sigma,V] = svd(X1,"econ");
S = U'*X2*V*diag(1./diag(Sigma)); %questa dovrebbe essere ATilde
[eV,D] = eig(S); %Eigenvectors and eigenvalues
mu = diag(D);
omega = log(mu)/(dt);
Phi = U*eV;
y0 = Phi\u; %pseudo-inverse initial conditions
u_modes = zeros(size(V,2),length(t));
for iter = 1:length(t)
    u_modes(:,iter) = y0.*exp(omega*t(iter));
end
u_dmd = Phi*u_modes;

%figure(3), waterfall(space,t',abs(u_dmd.'));
t_plot=inputdata(1,1:end);
figure(4), hold on, set(gca,'Fontsize',20), grid on, plot(t_plot,1000*pelts(1,:)','b-o','LineWidth',2),plot(t_plot,1000*abs(u_dmd(1,:)'),'c--','LineWidth',2), xlabel("Year"), ylabel("Population")
title('Snowshore Hare')
legend('Real Data','Exact DMD Approx.'), hold off;
figure(5), hold on, set(gca,'Fontsize',20), grid on, plot(t_plot,1000*pelts(2,:)','r-x','LineWidth',2),plot(t_plot,1000*abs(u_dmd(2,:)'),"Color",[0.4660 0.6740 0.1880],'LineStyle','--','LineWidth',2)
xlabel("Year"); ylabel("Population")
title('Lynx')
legend('Real Data','Exact DMD Approx.')
hold off;
% close all;
%% Time embedded
p = 12; contatore=0; r=3;
for h_idx = 1:p
    H(h_idx+contatore:h_idx*size(X,1),:) = X(:,1+contatore:end-p+contatore);
    contatore = contatore+1;
end

H1 = H(:,1:end-param);
H2 = H(:,param+1:end);
uh = H(:,1); %initial condition
[Uhtemp,Sigmahtemp,Vhtemp] = svd(H1,"econ");

figure
hold on
grid on
plot(sum(Sigmahtemp)./sum(sum(Sigmahtemp))*100)
hold off

figure
hold on
grid on
plot(Uhtemp(:,1:r))
hold off

figure
hold on
grid on
plot(Vhtemp(:,1:r))
hold off


Uh = Uhtemp(:,1:r);
Sigmah = Sigmahtemp(1:r,1:r);
Vh = Vhtemp(:,1:r);


Sh = Uh'*H2*Vh*diag(1./diag(Sigmah));
[eVh,Dh] = eig(Sh); %Eigenvectors and eigenvalues
muh = diag(Dh);
omegah = log(muh)/(dt);
Phih = Uh*eVh;
y0h = Phih\uh; %pseudo-inverse initial conditions
u_modesh = zeros(size(Vh,2),length(t));
for iter = 1:length(t)
    u_modesh(:,iter) = y0h.*exp(omegah*t(iter));
end
u_dmdh = Phih*u_modesh;

figure(6), hold on, grid on, plot(t_plot,pelts(1,:)'),plot(t_plot,abs(u_dmdh(1,:)')), hold off;
figure(7), hold on, grid on, plot(t_plot,pelts(2,:)'),plot(t_plot,abs(u_dmdh(2,:)')), hold off;