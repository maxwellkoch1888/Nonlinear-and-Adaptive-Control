clear all; close all; 
% error somewhere in code
y0 = [1;1;0];
[t,y] = ode45(@ode,[0,15], y0);

figure, plot(y(:,1),'linewidth',3)
grid on 

function ydot = ode(t,y) 

    x = y(1);
    xd = y(2); 
    sigma = y(3); 
    
    Kp = -1.5; 
    Kd = -2; 
    Ki = -1;

    d = 10; 
    
    % u = Kp*x + Kd*xd + Ki*sigma; 
    u = -xd - d*sign(x+xd);  
    
    xdot = xd; 
    xddot = -x + u + d*sin(t); 
    sigmadot = x;
    ydot = [xdot; xddot; sigmadot];

end 