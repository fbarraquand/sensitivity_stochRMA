function dydt=rma_holling(t,y)   % function handle for ODE

global r epsilon a b m K

dydt = zeros(3,1);   
dydt(1) = r*y(1)*(1.0 - y(1)/K) - a(1)*y(1)*y(2)/(1+a(1)*h(1)*y(1));
dydt(2) = epsilon*a(1)*y(1)*y(2)/(1+b(1)*y(1)) - m*y(2);
dydt(3) = 1;
