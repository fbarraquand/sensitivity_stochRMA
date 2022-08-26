%%% FB 31/03/2019 -- SDE version of RMA model to verify is structural sensitivity still works in a stochastic context
%%% This code has no functions (on purpose), because the integration algorithm uses brute force and hence we need to maximize speed
%%% 01/04/19 -- adds the tanh functional response

clc
clear all
close all

randn('seed',10)

%tic

% ------------------- Parameters --------------------------------%

r=1;
m=0.1;
epsilon=1;

a_H = 3.05;
b_H = 2.68;

a_I = 1;
b_I = 2;

a_T = 0.99;
b_T = 1.48;

x=linspace(0,10,1000)
figure,
subplot(121) 
hold on
plot(x,(a_H*x)./(1+b_H*x),'k-')
plot(x,a_I*(1-exp(-b_I*x)),'b-')
plot(x,a_T*tanh(b_T*x),'r-')
hold off
subplot(122) %% Something possibly smarter
C = 1
D = 0.5/log(2)
hold on
plot(x,(C*x)./(D+x),'k-')
plot(x,C*(1-exp(-log(2)*x/D)),'b-')
hold off
%% Ivlev saturates more slowly

% trigo functional response to implement

%---------------- Stochasticity ---------------------------------% 

%%% sigmaE = 0.5; % Intrinsic noise SD (I started with sigmaE=0.1 or 0.5)
sigmaE = 0.1
%sigmaE=0.05
%------------------- Numerical integration stuff ----------------% 
intervals = 20; %50 -> dt = 0.02, 20-> dt = 0.05
dt=1.0/intervals;  %-> NB decrease the timestep below theirs - they use RK4
T_max=2000; %50.0 or 8000
timevec=linspace(0,T_max,intervals*T_max+1);

%-------------------- Loop on K ---------------------------------%

%K=0.1:0.5:2.6;
K=0.1:0.005:1.5
length_stored=3000; %1500
%previously the last 1500 but 800/0.05 = 16 000, so this meant 75 time units and only 3.75 periods on average around the bifurcation point. Not much considering the period can increase. 
% Hence we now pick 1 every 10 time points (see below) of 30000 points which means 1500 time units (the 500 first time units are a warm-up to remove transients). 
% This helps to compute quantiles. 
nmatrix=zeros(length(K),length_stored);
pmatrix=zeros(length(K),length_stored);
n2matrix=zeros(length(K),length_stored);
p2matrix=zeros(length(K),length_stored);

nmatrix_I=zeros(length(K),length_stored);
pmatrix_I=zeros(length(K),length_stored);
n2matrix_I=zeros(length(K),length_stored);
p2matrix_I=zeros(length(K),length_stored);

nmatrix_T=zeros(length(K),length_stored);
pmatrix_T=zeros(length(K),length_stored);
n2matrix_T=zeros(length(K),length_stored);
p2matrix_T=zeros(length(K),length_stored);

for (kk=1:length(K)),

disp(kk)

%%% That's what would be logical but here epsilon=1 goes against logic
%n(1) = 0.1;
%p(1) = 0.05;
%n2(1) = 0.1;
%p2(1) = 0.05;

n(1) = 0.5;
p(1) = 0.5;
n2(1) = 0.5;
p2(1) = 0.5;

n_I(1) = 0.5;
p_I(1) = 0.5;
n2_I(1) = 0.5;
p2_I(1) = 0.5;

n_T(1) = 0.5;
p_T(1) = 0.5;
n2_T(1) = 0.5;
p2_T(1) = 0.5;


%-------------------- Loop on time ------------------------------%
ki=1;
t=0.0;
while(t<(T_max-dt))

PreyGrowth = r*n(ki)*(1-n(ki)/K(kk));
Predation = (a_H*n(ki)*p(ki))/(b_H*n(ki)+1);
PredGrowth = epsilon*Predation-m*p(ki);

PreyGrowth2 = r*n2(ki)*(1-n2(ki)/K(kk));
Predation2 = (a_H*n2(ki)*p2(ki))/(b_H*n2(ki)+1);
PredGrowth2 = epsilon*Predation2-m*p2(ki);

PreyGrowth_I = r*n_I(ki)*(1-n_I(ki)/K(kk));
Predation_I = a_I*(1-exp(-b_I*n_I(ki)))*p_I(ki);
PredGrowth_I = epsilon*Predation_I-m*p_I(ki);

PreyGrowth2_I = r*n2_I(ki)*(1-n2_I(ki)/K(kk));
Predation2_I = a_I*(1-exp(-b_I*n2_I(ki)))*p2_I(ki);
PredGrowth2_I = epsilon*Predation2_I-m*p2_I(ki);

PreyGrowth_T = r*n_T(ki)*(1-n_T(ki)/K(kk));
Predation_T = a_T*tanh(b_T*n_T(ki))*p_T(ki);
PredGrowth_T = epsilon*Predation_T-m*p_T(ki);

PreyGrowth2_T = r*n2_T(ki)*(1-n2_T(ki)/K(kk));
Predation2_T =  a_T*tanh(b_T*n2_T(ki))*p2_T(ki);
PredGrowth2_T = epsilon*Predation2_T-m*p2_T(ki);


%----------------------------------------------------%
%%% This is Euler-Maruyama (simplest) scheme 
%%% Multiplication by abundance for environmental stochasticity (instead of sqrt(abundance) for dem. stoch. see S. Engen papers)
% Env. stochasticity
StochCompPrey = sigmaE*n(ki)*sqrt(dt)*randn; % Standard deviation of the Wiener process increment is sqrt(dt)
StochCompPred = sigmaE*p(ki)*sqrt(dt)*randn; % We use the same sigmaE for parsimony, though there's no other reason
StochCompPrey_I = sigmaE*n_I(ki)*sqrt(dt)*randn; % Standard deviation of the Wiener process increment is sqrt(dt)
StochCompPred_I = sigmaE*p_I(ki)*sqrt(dt)*randn;
StochCompPrey_T = sigmaE*n_T(ki)*sqrt(dt)*randn; % Standard deviation of the Wiener process increment is sqrt(dt)
StochCompPred_T = sigmaE*p_T(ki)*sqrt(dt)*randn;

%----------------------------------------------------%
n(ki+1) = n(ki) + dt*(PreyGrowth - Predation) + StochCompPrey;%+ DemogCompPrey; % Only env. stoch. on the prey
p(ki+1) = p(ki) + dt*PredGrowth + StochCompPred ; % + DemogCompPred Put only dem. stoch. on the predator, ideally

% deterministic equivalent
n2(ki+1) = n2(ki) + dt*(PreyGrowth2 - Predation2); % 
p2(ki+1) = p2(ki) + dt*PredGrowth2; %

n_I(ki+1) = n_I(ki) + dt*(PreyGrowth_I - Predation_I) + StochCompPrey_I;%+ DemogCompPrey; % Only env. stoch. on the prey
p_I(ki+1) = p_I(ki) + dt*PredGrowth_I + StochCompPred_I ; % + DemogCompPred Put only dem. stoch. on the predator, ideally

% deterministic equivalent
n2_I(ki+1) = n2_I(ki) + dt*(PreyGrowth2_I - Predation2_I); % 
p2_I(ki+1) = p2_I(ki) + dt*PredGrowth2_I; %

n_T(ki+1) = n_T(ki) + dt*(PreyGrowth_T - Predation_T) + StochCompPrey_T;%+ DemogCompPrey; % Only env. stoch. on the prey
p_T(ki+1) = p_T(ki) + dt*PredGrowth_T + StochCompPred_T ; % + DemogCompPred Put only dem. stoch. on the predator, ideally

% deterministic equivalent
n2_T(ki+1) = n2_T(ki) + dt*(PreyGrowth2_T - Predation2_T); % 
p2_T(ki+1) = p2_T(ki) + dt*PredGrowth2_T; %


ki=ki+1;
t=t+dt;
end
%-----------------------end of time loop --------------------------%

%% plot the time series
%if (0==1), %should forbid plotting
if (K(kk)==0.1 | K(kk)==0.5 | K(kk)==1), 
set(gca,'FontSize',7,'FontName','Arial')

figure,
subplot(321)
plot(timevec,n,timevec,n2,'k--','LineWidth',1)
ylabel('\fontsize{7}Prey')
xlabel('\fontsize{7}Time')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
text(0,1.25,'(a)','Units','normalized','Fontsize',8)
axis tight

subplot(322)
plot(timevec,p,timevec,p2,'k--','LineWidth',1)
ylabel('\fontsize{7}Predator')
legend('\fontsize{7}SDE', '\fontsize{7}ODE','Location','northeast')% 'northeastoutside' shrinks the x-axis and best is not really best
axis tight
xlabel('\fontsize{7}Time')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
title('\fontsize{7}Holling')
text(0,1.25,'(b)','Units','normalized','Fontsize',8)

subplot(323)
plot(timevec,n_I,timevec,n2_I,'k--','LineWidth',1)
ylabel('\fontsize{7}Prey')
xlabel('\fontsize{7}Time')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
text(0,1.25,'(c)','Units','normalized','Fontsize',8)
axis tight

subplot(324)
plot(timevec,p_I,timevec,p2_I,'k--','LineWidth',1)
ylabel('\fontsize{7}Predator')
%legend('SDE', 'ODE')
axis tight
xlabel('\fontsize{7}Time')
title('\fontsize{7}Ivlev')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
text(0,1.25,'(d)','Units','normalized','Fontsize',8)

subplot(325)
plot(timevec,n_T,timevec,n2_T,'k--','LineWidth',1)
ylabel('\fontsize{7}Prey')
xlabel('\fontsize{7}Time')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
text(0,1.25,'(e)','Units','normalized','Fontsize',8)
axis tight

subplot(326)
plot(timevec,p_T,timevec,p2_T,'k--','LineWidth',1)
ylabel('\fontsize{7}Predator')
%legend('SDE', 'ODE')
axis tight
xlabel('\fontsize{7}Time')
title('\fontsize{7}tanh')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
text(0,1.25,'(f)','Units','normalized','Fontsize',8)

print(figure(gcf),'-dpng','-r600',['RMA_Traj_K',num2str(K(kk))]);
end


%% store the end of it - now every 10 points
nmatrix(kk,:)=n((end-10*length_stored+1):10:end);
pmatrix(kk,:)=p((end-10*length_stored+1):10:end);
n2matrix(kk,:)=n2((end-10*length_stored+1):10:end);
p2matrix(kk,:)=p2((end-10*length_stored+1):10:end);

nmatrix_I(kk,:)=n_I((end-10*length_stored+1):10:end);
pmatrix_I(kk,:)=p_I((end-10*length_stored+1):10:end);
n2matrix_I(kk,:)=n2_I((end-10*length_stored+1):10:end);
p2matrix_I(kk,:)=p2_I((end-10*length_stored+1):10:end);

nmatrix_T(kk,:)=n_T((end-10*length_stored+1):10:end);
pmatrix_T(kk,:)=p_T((end-10*length_stored+1):10:end);
n2matrix_T(kk,:)=n2_T((end-10*length_stored+1):10:end);
p2matrix_T(kk,:)=p2_T((end-10*length_stored+1):10:end);
end

%-----------------------end of K loop -----------------------------%
figure,
set(gca,'FontSize',7,'FontName','Arial'); %otherwise Helvetica

subplot(321)
plot(K,log(nmatrix),'k.','MarkerSize',2) 
hold on
plot(K,log(n2matrix),'r.','MarkerSize',2)
xlabel('')
ylabel('\fontsize{7}log(prey density)')
text(0,1.25,'(a)','Units','normalized','Fontsize',7)
axis([min(K) 1.25 -45 5])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
title('\fontsize{7}Holling')
hold off


subplot(322)
plot(K,nmatrix,'k.','MarkerSize',2)
%plot(K,nmatrix,'ko','MarkerFaceColor','k','MarkerSize',2)
hold on
plot(K,quantile(nmatrix.',0.975),'Color',[0.7 0.7 0.7]) 
plot(K,n2matrix,'r.','MarkerSize',2)
xlabel('')
ylabel('\fontsize{7}prey density')
text(0,1.25,'(b)','Units','normalized','Fontsize',7)
axis([min(K) 1.25 0 2.5])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
title('\fontsize{7}Holling')
hold off


subplot(323)
plot(K,log(nmatrix_I),'k.','MarkerSize',2)
hold on
plot(K,log(n2matrix_I),'r.','MarkerSize',2)
xlabel('')
ylabel('\fontsize{7}log(prey density)')
text(0,1.25,'(c)','Units','normalized','Fontsize',7)
axis([min(K) 1.5 -45 5])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
title('\fontsize{7}Ivlev')
hold off

subplot(324)
plot(K,nmatrix_I,'k.','MarkerSize',2)
hold on
plot(K,quantile(nmatrix_I.',0.975),'Color',[0.7 0.7 0.7]) 
plot(K,n2matrix_I,'r.','MarkerSize',2)
xlabel('')
ylabel('\fontsize{7}prey density')
text(0,1.25,'(d)','Units','normalized','Fontsize',7)
axis([min(K) 1.5 0 2.5])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
title('\fontsize{7}Ivlev')

subplot(325)
plot(K,log(nmatrix_T),'k.','MarkerSize',2)
hold on
plot(K,log(n2matrix_T),'r.','MarkerSize',2)
xlabel('\fontsize{7}K')
ylabel('\fontsize{7}log(prey density)')
text(0,1.25,'(e)','Units','normalized','Fontsize',7)
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
axis([min(K) 1.5 -45 5])
title('\fontsize{7}tanh')
hold off

subplot(326)
plot(K,nmatrix_T,'k.','MarkerSize',2)
hold on
plot(K,quantile(nmatrix_T.',0.975),'Color',[0.7 0.7 0.7]) 
plot(K,n2matrix_T,'r.','MarkerSize',2)
xlabel('\fontsize{7}K')
ylabel('\fontsize{7}prey density')
text(0,1.25,'(f)','Units','normalized','Fontsize',7)
axis([min(K) 1.5 0 2.5])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',7);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',7);
title('\fontsize{7}tanh')
hold off

%toc
%
print(figure(gcf),'-dpng','-r300','RMA_Kloop_sigma01');
