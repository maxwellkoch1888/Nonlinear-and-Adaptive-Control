clear all; close all; clc 

% Solve the diff eq yddot = cos(2y)
y0 = [0,1.1];

[t,y] = ode45(@ode,[0,15],y0);
figure, plot(t,y(:,1), 'linewidth',3), grid on
xlabel('t'), ylabel('y1')

function ydot = ode(t,y)
    y1 = y(1);
    y2 = y(2); 
    y1dot = y2;
    y2dot = cos(2*y1);
    ydot = [y1dot; y2dot];
end