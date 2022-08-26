function dydt=rma_tanh(t,y)   % function handle for ODE

global r epsilon a b m K

dydt = zeros(3,1);   
dydt(1) = r*y(1)*(1.0 - y(1)/K) - a(3)*tanh(b(3)*y(1))*y(2);
dydt(2) = epsilon*a(3)*tanh(b(3)*y(1))*y(2) - m*y(2);
dydt(3) = 1;
