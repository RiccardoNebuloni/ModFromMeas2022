clear; close all; clc;

td = [0:2:29*2];

hare = [20 20 52 83 64 68 83 12 36 150 110 60 7 10 70 100 92 70 10 11 137 137 18 22 52 83 18 10 9 65] ;
lynx = [32 50 12 10 13 36 15 12 6 6 65 70 40 9 20 34 45 40 15 15 60 80 26 18 37 50 35 12 12 25];


%p = [20 32  0.6156    0.0301    0.5706    0.0114];
 %p= [initial population(1) initial population(2) parameters(1),parameters(2),parameters(3), parameters(4)]
 p = [hare(1) lynx(1) 0.6 0.0301 0.5706 0.0114];


% fminsearch routine 

fun = @leastcomp;
[p,fval,exitflag] = fminsearch(@leastcomp,p,[],td,hare,lynx);
p
fval

[t,y] = ode45(@lotvol,td,[p(1),p(2)],[],p(3),p(4),p(5),p(6));

% figure
% plot(t,y)

%% confronto
figure(1)
hold on
grid on
title 'Shoe'
plot(t,hare',LineWidth=2);
plot(t,y(:,1),LineWidth=2);
legend('Real','Fitting')
hold off

figure(2)
hold on
grid on
title 'Lynx'
plot(t,lynx',LineWidth=2);
plot(t,y(:,2),LineWidth=2);
legend('Real','Fitting')
hold off
