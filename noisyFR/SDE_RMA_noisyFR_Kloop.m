%%% FB 31/03/2019 -- SDE version of RMA model to verify is structural sensitivity still works in a stochastic context


clc
clear all
close all

randn('seed',10)

% ------------------- Parameters --------------------------------%

r=1;
m=0.1;
epsilon=1;

a_H = 3.05;
b_H = 2.68;

a_I = 1;
b_I = 2;

gamma = 0.5 % autocorrelation FR

x=linspace(0,10,1000)
figure,
subplot(121) 
hold on
plot(x,(a_H*x)./(1+b_H*x),'k-')
plot(x,a_I*(1-exp(-b_I*x)),'b-')
hold off
subplot(122) %% Something smarter
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
%sigmaE = 0.25
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

noiseFR = 0;

%-------------------- Loop on time ------------------------------%
ki=1;
t=0.0;
while(t<(T_max-dt))

noiseFR = -gamma*noiseFR*dt + 0.25*sqrt(dt)*randn; %% Noise on the FR, on the log(half saturation) -> thus multiplicative on the half-sat scale

PreyGrowth = r*n(ki)*(1-n(ki)/K(kk));
Predation = (C*n(ki)*p(ki))/(n(ki)+D*exp(noiseFR));
PredGrowth = epsilon*Predation-m*p(ki);

PreyGrowth2 = r*n2(ki)*(1-n2(ki)/K(kk));
Predation2 = (C*n2(ki)*p2(ki))/(n2(ki)+D);
PredGrowth2 = epsilon*Predation2-m*p2(ki);

PreyGrowth_I = r*n_I(ki)*(1-n_I(ki)/K(kk));
Predation_I = C*(1-exp(-log(2)*n_I(ki)/(D*exp(noiseFR) ) ))*p_I(ki);
PredGrowth_I = epsilon*Predation_I-m*p_I(ki);

PreyGrowth2_I = r*n2_I(ki)*(1-n2_I(ki)/K(kk));
Predation2_I = C*(1-exp(-log(2)*n2_I(ki)/(D)))*p2_I(ki);
PredGrowth2_I = epsilon*Predation2_I-m*p2_I(ki);

%----------------------------------------------------%
%%% This is Euler-Maruyama (simplest) scheme 
%%% Multiplication by abundance for environmental stochasticity (instead of sqrt(abundance) for dem. stoch. see S. Engen papers)
% Env. stochasticity
StochCompPrey = sigmaE*n(ki)*sqrt(dt)*randn; % Standard deviation of the Wiener process increment is sqrt(dt)
StochCompPred = sigmaE*p(ki)*sqrt(dt)*randn; % We use the same sigmaE for parsimony, though there's no other reason
StochCompPrey_I = sigmaE*n_I(ki)*sqrt(dt)*randn; % Standard deviation of the Wiener process increment is sqrt(dt)
StochCompPred_I = sigmaE*p_I(ki)*sqrt(dt)*randn;
% Dem. stochasticity 
%DemogCompPrey = sqrt(abs(PreyGrowth) + PredGen + PredSpe)*sqrt(dt)*randn; % Demographic stochasticity - maybe problematic on density rather than abundance
%%DemogCompPred = sqrt(s*p(ki))*sqrt(dt)*randn; %%% sqrt(abs(PredGrowth)) a bit dirty - better to separate birth death rates, and sum them... but there we mostly have birth
%%% Some problems with ddem. stochasticity in the predator - sqrt(p) with p<1 yields large fluctuations, problems of scaling...
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


ki=ki+1;
t=t+dt;
end
%-----------------------end of time loop --------------------------%

%% plot the time series
if (K(kk)==0.1 | K(kk)==0.5 | K(kk)==1), 
figure,
set(gca,'FontSize',8)

subplot(221)
plot(timevec(1:5000),n(1:5000),timevec(1:5000),n2(1:5000),'k--','LineWidth',1)
ylabel('\fontsize{8}Prey')
xlabel('\fontsize{8}Time')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
text(0,1.1,'(a)','Units','normalized','Fontsize',8)
axis tight

subplot(222)
plot(timevec(1:5000),p(1:5000),timevec(1:5000),p2(1:5000),'k--','LineWidth',1)
ylabel('\fontsize{8}Predator')
legend('\fontsize{8}SDE', '\fontsize{8}ODE','Location','northeast')% 'northeastoutside' shrinks the x-axis and best is not really best
axis tight
xlabel('\fontsize{8}Time')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
title('\fontsize{8}Holling')
text(0,1.1,'(b)','Units','normalized','Fontsize',8)

subplot(223)
plot(timevec(1:5000),n_I(1:5000),timevec(1:5000),n2_I(1:5000),'k--','LineWidth',1)
ylabel('\fontsize{8}Prey')
xlabel('\fontsize{8}Time')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
text(0,1.1,'(c)','Units','normalized','Fontsize',8)
axis tight

subplot(224)
plot(timevec(1:5000),p_I(1:5000),timevec(1:5000),p2_I(1:5000),'k--','LineWidth',1)
ylabel('\fontsize{8}Predator')
%legend('SDE', 'ODE')
axis tight
xlabel('\fontsize{8}Time')
title('\fontsize{8}Ivlev')
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
text(0,1.1,'(d)','Units','normalized','Fontsize',8)


print(figure(gcf),'-dpng','-r300',['RMA_noisyFR_Traj_K',num2str(K(kk))]);
end

nmatrix(kk,:)=n((end-10*length_stored+1):10:end);
pmatrix(kk,:)=p((end-10*length_stored+1):10:end);
n2matrix(kk,:)=n2((end-10*length_stored+1):10:end);
p2matrix(kk,:)=p2((end-10*length_stored+1):10:end);

nmatrix_I(kk,:)=n_I((end-10*length_stored+1):10:end);
pmatrix_I(kk,:)=p_I((end-10*length_stored+1):10:end);
n2matrix_I(kk,:)=n2_I((end-10*length_stored+1):10:end);
p2matrix_I(kk,:)=p2_I((end-10*length_stored+1):10:end);
end
%-----------------------end of K loop -----------------------------%

figure,
set(gca,'FontSize',8)

subplot(221)
plot(K,log(nmatrix),'k.','MarkerSize',2)
hold on
plot(K,log(n2matrix),'r.','MarkerSize',2)
xlabel('K')
ylabel('\fontsize{8}log(prey density)')
text(0,1.1,'(a)','Units','normalized','Fontsize',8)
axis([0 1.5 -45 5])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
title('\fontsize{8}Holling')


subplot(222)
plot(K,nmatrix,'k.','MarkerSize',2)
hold on
plot(K,n2matrix,'r.','MarkerSize',2)
xlabel('\fontsize{8}K')
ylabel('\fontsize{8}prey density')
text(0,1.1,'(b)','Units','normalized','Fontsize',8)
axis([0 1.5 0 3])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
title('\fontsize{8}Holling')


subplot(223)
plot(K,log(nmatrix_I),'k.','MarkerSize',2)
hold on
plot(K,log(n2matrix_I),'r.','MarkerSize',2)
xlabel('\fontsize{8}K')
ylabel('\fontsize{8}log(prey density)')
text(0,1.1,'(d)','Units','normalized','Fontsize',8)
axis([0 1.5 -45 5])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
title('\fontsize{8}Ivlev')

subplot(224)
plot(K,nmatrix_I,'k.','MarkerSize',2)
hold on
plot(K,n2matrix_I,'r.','MarkerSize',2)
xlabel('\fontsize{8}K')
ylabel('\fontsize{8}prey density')
text(0,1.1,'(d)','Units','normalized','Fontsize',8)
axis([0 1.5 0 3])
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
title('\fontsize{8}Ivlev')

print(figure(gcf),'-dpng','-r300','RMA_noisyFR_Kloop');

