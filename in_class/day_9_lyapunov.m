clear all; close all; 

eq = 20*randn(3,1);

sol = fsolve(@fun, eq)

function F = fun(x)
x1 = x(1); 
x2 = x(2); 
x3 = x(3);

F = [-x1 + x2; -1/10*x1^3 - 10*sin(x3); x2 - x3];
end