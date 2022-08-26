%% test xstar


clc
clear all
close all

global a b r

r = 1;
a = [3.05 1 0.99];
b = [2.68 2 1.48];
m = 0.1;
K = 0.1:0.01:1.5;

func_resp_inverse(m,1)
func_resp_inverse(m,2)
func_resp_inverse(m,3)

function y = func_resp_inverse(x,model_type)
global a b
switch model_type 
  case 1
   y = a(1)*b(1)/(x-a(1));
  case 2
   y = (-1/b(2))*log((a(2)-x)/a(2))
  case 3
   y = (1/b(3))*atanh(x/a(3));
end
end



