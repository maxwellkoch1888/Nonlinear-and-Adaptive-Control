function [u, zdot] = pi_controller(t,state,x_eq,u_eq,model)
% omega_dot = f(x) + G(x)*controls

% Pull out actual rates
p = state(4); 
q = state(5); 
r = state(6); 
omega = state(4:6);

% Define reference signal
omega_ref = pi/180 * [5*sin(2*t); -2*sin(t); sin(t)]; 
omega_ref_dot = pi/180 * [10*cos(2*t); -2*cos(t); cos(t)]; 
error = omega - omega_ref;

% Define error
int_error = state(13:15); % integral error

% Define desired angular acceleration with proportional controller
Kp = diag([4,2,4]);
Ki = diag([1,1,1]); 
omegadot_des = -Kp*error - Ki*int_error + omega_ref_dot; 

% Current model output
[~, ~, hmat, Imat, ~, M_b] = model(state,u_eq); 

% Pull out inertia terms, calculate inertia term
Ixx = Imat(1,1);  Iyy = Imat(2,2);  Izz = Imat(3,3);
Ixy = -Imat(1,2); Ixz = -Imat(1,3); Iyz = -Imat(2,3);

inertia_effects = [(Iyy-Izz)*q*r + Iyz*(q^2-r^2) + Ixz*p*q - Ixy*p*r; 
                   (Izz-Ixx)*p*r + Ixz*(r^2-p^2) + Ixy*q*r - Iyz*p*q;
                   (Ixx-Iyy)*p*q + Ixy*(p^2-q^2) + Iyz*p*r - Ixz*q*r];

% f(x) 
f = Imat \ (M_b + hmat*omega + inertia_effects);

% Calculate G(x) numerically
eps = 1e-6;
G = zeros(3,3);

for i = 1:3
    du = zeros(4,1);
    du(i) = eps;

    [~,~,~,~,~,M_pert] = model(state,u_eq + du);
    f_pert = Imat \ (M_pert + hmat*omega + inertia_effects);

    G(:,i) = (f_pert - f)/eps;
end

% Dynamic inversion
delta = pinv(G)*(omegadot_des - f);

% Build full control vector
u = u_eq;
u(1:3) = u_eq(1:3) + delta;
u(4) = 0.282388376195832;
zdot = error;
end

