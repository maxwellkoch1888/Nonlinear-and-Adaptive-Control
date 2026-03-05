clear all; close all; 

% Trim condition and model 
x_eq = 1e4*[0.063502542642571; 0; 0.003196923248873; 0; 0; 0; 0; 0.000005030076610; 0; 0; 0; -1.500000000000000];
u_eq = [0; -0.001389023473370; 0; 0.282388376195832];
model = @base_aero;

% Define controller type: 1 for p, 2 for reference tracking
ctrl = 1;

% Initial state with perturbation
x0 = x_eq;
x0(4) = x0(4) - 5*pi/180;
x0(5) = x0(5) + 1*pi/180;
x0(6) = x0(6) + 3*pi/180;
if ctrl == 2 
    x0 = [x0; 0; 0; 0];   % add integral states
end 

% Run simulation
ops = odeset('AbsTol',1e-10,'RelTol',1e-10, 'OutputFcn', @odeplot);
[t,x] = ode45(@(t,state) ode_0th(t,state,model,x_eq,u_eq,ctrl), ...
              [0 5], x0, ops);

if ctrl == 2
    omega_ref_plot = pi/180 * [
        5*sin(2*t), ...
       -2*sin(t), ...
        sin(t)
    ];
end 

% Plot results
figure

subplot(2,2,1)
plot(t, x(:,1), 'LineWidth', 2); hold on
plot(t, x(:,2), 'LineWidth', 2);
plot(t, x(:,3), 'LineWidth', 2);
grid on
title('Velocities')
ylabel('Velocity [ft/sec]')
xlabel('Time [sec]')
legend('u','v','w')

subplot(2,2,2)
plot(t, x(:,4)*180/pi, 'LineWidth', 2); hold on
plot(t, x(:,5)*180/pi, 'LineWidth', 2);
plot(t, x(:,6)*180/pi, 'LineWidth', 2);
if ctrl == 2
    plot(t, omega_ref_plot(:,1)*180/pi, '--', 'LineWidth', 2)
    plot(t, omega_ref_plot(:,2)*180/pi, '--', 'LineWidth', 2)
    plot(t, omega_ref_plot(:,3)*180/pi, '--', 'LineWidth', 2)
end 
grid on
title('Angular Rates')
ylabel('Rate [deg]')
xlabel('Time [sec]')
if ctrl == 2 
    legend('p','q','r', 'p_{ref}', 'q_{ref}', 'r_{ref}')
else 
    legend('p','q','r')
end 

subplot(2,2,3)
plot(t, x(:,7)*180/pi, 'LineWidth', 2); hold on
plot(t, x(:,8)*180/pi, 'LineWidth', 2);
plot(t, x(:,9)*180/pi, 'LineWidth', 2);
grid on
title('Euler Angles')
ylabel('Angle [deg]')
xlabel('Time [sec]')
legend('\phi','\theta','\psi')

subplot(2,2,4)
plot(t, x(:,10), 'LineWidth', 2); hold on
plot(t, x(:,11), 'LineWidth', 2);
plot(t, x(:,12), 'LineWidth', 2);
grid on
title('Position')
ylabel('Position [ft]')
xlabel('Time [sec]')
legend('x_f','y_f','z_f')