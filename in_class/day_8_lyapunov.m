clear all; close all; 

x0 = [1;1];
ops = ['AbsTol',1e-12,'RelTol',1e-12];
[t,x] = ode45(@ode, [0,100], x0);

figure, plot(x(:,1), x(:,2), 'linewidth', 2), grid on 

V = x(:,1).^4/4 + x(:,2).^2/2;
% figure, plot(t,V)

function xdot = ode(t,x)
    
    x1 = x(1); 
    x2 = x(2); 
    u = tanh(100*x2);
    x1dot = x2; 
    x2dot = -x1^3-x2^2*u; 
    xdot = [x1dot; x2dot];

end 