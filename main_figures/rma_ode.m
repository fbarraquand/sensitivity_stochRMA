%% Numerical integration RMA model


clc
clear all
close all

global r a b epsilon K m

r=1;
epsilon=1;
K=5;
m=0.1;

a = [3.05 1 0.99];
b = [2.68 2 1.48];

y0=[0.5 0.5 0];                    % initial values of state variables vector: plant - prey - predator - time 
tstart = 0;
tstop = 1150; %30, 250
tspan=[tstart tstop];                        % timespan for the numerical solution

options= odeset('Reltol',1e-4,'NonNegative',[1 2]);
[tout,yout] = ode45(@rma_tanh , tspan , y0,options); 

% Plotting
subplot(211)
plot(tout,yout(:,1))
subplot(212)
plot(tout,yout(:,2))


%%% Just below

K=3;
m=0.1;

a = [3.05 1 0.99];
b = [2.68 2 1.48];

y0=[0.5 0.5 0];                    % initial values of state variables vector: plant - prey - predator - time 
tstart = 0;
tstop = 1150; %30, 250
tspan=[tstart tstop];                        % timespan for the numerical solution

options= odeset('Reltol',1e-4,'NonNegative',[1 2]);
[tout,yout] = ode45(@rma_tanh , tspan , y0,options); 

% Plotting
subplot(211)
plot(tout,yout(:,1))
subplot(212)
plot(tout,yout(:,2))

