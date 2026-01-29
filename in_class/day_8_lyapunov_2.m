clear all; close all; 

x0 = 1;
ops = ['AbsTol',1e-12,'RelTol',1e-12];
[t,x] = ode45(@ode, [0,10], x0);
u = -x./sin(x); 
figure, plot(t,x), grid on 
figure, plot(t,u)

% V = x(:,1).^4/4;
% figure, plot(t,V)

function xdot = ode(t,x)
    

    u = -x/sin(x);
    xdot = sin(x)*u; 

end 