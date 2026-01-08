clear all; close all;

x0 = [3;3];
[t,x] = ode45(@ode,[0,25],x0);
figure, plot(x(:,1),x(:,2), 'linewidth',3)
xlabel('x1'), ylabel('x2')

A = [0, 1; -1,1];
eig(A)

function xdot = ode(t,x)
    x1 = x(1);
    x2 = x(2);
    x1dot = x2;
    x2dot = -x1 + (1-x1^2)*x2;
    xdot = [x1dot;x2dot];
end