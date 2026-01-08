clear all; close all; clc

% Data 
L = 1;
g = 10;
m = 1;

[t,x] = ode45(@ode,[0,5], [45*pi/180;0],[],L,g);

figure, plot(t,x(:,1))

function xdot = ode(t,x,L,g)

x1 = x(1);
x2 = x(2);
u = -10 * x2; 

x1dot = x2;
x2dot = -g/L * sin(x1) + u;
xdot = [x1dot; x2dot];

end
