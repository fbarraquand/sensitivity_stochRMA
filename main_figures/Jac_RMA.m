%%% Eigenvalues properties of the RMA model -- always a stable focus? 
%%% FB 01/04/2019 -- currently just over a gradient of K

clc
clear all
close all

global a b r

r = 1
a = [3.05 1 0.99];
b = [2.68 2 1.48];
m = 0.1;
K = 0.01:0.01:1.5;

rp=zeros(length(K),3);
ip=zeros(length(K),3);

%% Loop on K
for kk=1:length(K),

 for model_type=1:3

xstar = func_resp_inverse(m,model_type);
ystar = g(xstar,K(kk))/func_resp(xstar,model_type);
if(imag(xstar)~=0), disp(['Equilibrium amiss',kk]); end

Jac = [gprime(xstar,K(kk))-ystar*func_resp_prime(xstar,model_type),  -func_resp(xstar,model_type); ystar*func_resp_prime(xstar,model_type), func_resp(xstar,model_type)-m];

real_part = real(eig(Jac));
[y,I] = max(abs(real_part));

rp(kk,model_type) = real_part(I);
ip(kk,model_type) = max(imag(eig(Jac)));

end

end

set(gca,'FontSize',8,'FontName','Arial'); %otherwise Helvetica

subplot(311) %plot eigenvalue real part
hold on 
plot(K,rp(:,1),'b-','LineWidth',2)
axis([0 1.5 -0.05 0.1])
pbaspect([3 1 1])
%ylabel('Re(\lambda)')
xpos = -0.25 % https://stackoverflow.com/questions/10634923/align-the-ylabel-in-subplots
yl=ylabel('Re(\lambda)');
pos=get(yl,'Pos')
set(yl,'Pos',[xpos pos(2) pos(3)])
title('\fontsize{8}Real part (zoom)')
text(0,1.2,'(a)','Units','normalized','FontSize',8)
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
plot(K,rp(:,2),'k-','LineWidth',2)
plot(K,rp(:,3),'r-','LineWidth',2)
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = '--';
hold off
%legend('\fontsize{10}Holling','\fontsize{10}Ivlev','\fontsize{10}Tanh','Location','northwest')

subplot(312) %plot eigenvalue real part
hold on 
plot(K,rp(:,1),'b-','LineWidth',2)
axis([0 1.5 -3 0.1])
pbaspect([3 1 1])
yl=ylabel('Re(\lambda)');
pos=get(yl,'Pos')
set(yl,'Pos',[xpos pos(2) pos(3)])
title('\fontsize{8}Real part')
text(0,1.2,'(b)','Units','normalized','FontSize',8)
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
plot(K,rp(:,2),'k-','LineWidth',2)
plot(K,rp(:,3),'r-','LineWidth',2)
%xline(0.2,'k--')
h=line([0.15 0.15], [-3 0.1]);
set(h,'Color','k','LineStyle','--')
hold off
%legend('\fontsize{10}Holling','\fontsize{10}Ivlev','\fontsize{10}Tanh')


subplot(313) % plot eigenvalue imaginary part [do I divide it by 2*pi to give a frequency?]
hold on 
plot(K,ip(:,1),'b-','LineWidth',2)
axis([0 1.5 0 0.31])
pbaspect([3 1 1])
xlabel('K')
yl=ylabel('Im(\lambda)');
pos=get(yl,'Pos')
set(yl,'Pos',[xpos pos(2) pos(3)])
title('\fontsize{8}Imaginary part')
text(0,1.2,'(c)','Units','normalized','FontSize',8)
xAX = get(gca,'XAxis');
set(xAX,'FontSize',8);
yAX = get(gca,'YAxis');
set(yAX,'FontSize',8);
plot(K,ip(:,2),'k-','LineWidth',2)
plot(K,ip(:,3),'r-','LineWidth',2)
h=line([0.15 0.15], [-3 0.35]);
set(h,'Color','k','LineStyle','--')
legend('\fontsize{8}Holling','\fontsize{8}Ivlev','\fontsize{8}tanh')
hold off

print(figure(gcf),'-dpng','-r600','EigenvaluesRMA.png')

%% This shows that 


%% Functions definitions

function y = g(x,K0)
global r
y = r*x - r*x^2/K0
end


function y = gprime(x,K0)
global r
y = r - 2*r*x/K0
end

function y = func_resp(x,model_type)
global a b
switch model_type 
  case 1
   y = a(1)*x/(1+b(1)*x);
  case 2
   y = a(2)*( 1-exp(-b(2)*x) );
  case 3
   y = a(3)*tanh(b(3)*x);
end
end

function y = func_resp_prime(x,model_type)
global a b
switch model_type 
  case 1
   y = a(1)/((1+b(1)*x)^2);
  case 2
   y = a(2)*b(2)*exp(-b(2)*x);
  case 3
   y = a(3)*b(3)*(1-(tanh(b(3)*x)^2));
end
end

function y = func_resp_inverse(x,model_type)
global a b
switch model_type 
  case 1
   y = x/(a(1)-b(1)*x);
  case 2
   y = (-1/b(2))*log((a(2)-x)/a(2))
  case 3
   y = (1/b(3))*atanh(x/a(3));
end
end




