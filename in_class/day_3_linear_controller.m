clear all; close all;

x0 = 10; 
sigma0 = 0;
y0 = [x0; sigma0];
[t,y] = ode45(@ode, [0,10], y0); 
figure, plot(t,y(:,1),'linewidth',3)
grid on 
xlabel('t'), ylabel('x')

function ydot = ode(t,y) 

    x = y(1);
    sigma = y(2);

    a = 7; 
    b = 1; 
    d = 1;

    Kp = -10; 
    Ki = -10;
    u = Kp*x + Ki*sigma;
    
    % Integral control can catch sinusoidal disturbance if frequency/ magnitude relation is small enough
    % xdot = a*x + b*u + d*sin(0.1*t);
    % Add a frictional distrurbance, nonlinear, can still catch if small enough
    % xdot = a*x + b*u + d*x^2;
    sigma_dot = x; 

    ydot = [xdot; sigma_dot]; 

end 