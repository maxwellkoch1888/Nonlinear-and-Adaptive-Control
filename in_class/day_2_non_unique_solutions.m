clear all; close all; clc 

%Demonstrate non unique solutions
% x0 = 0;
x0 = 1e-8;

[t,x] = ode45(@ode,[0,1],x0);
figure, plot(t,x)

function xdot = ode(t,x)
    xdot = x^(1/3);
end