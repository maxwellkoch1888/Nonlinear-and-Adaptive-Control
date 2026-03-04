function [state_dot,u] = ode_0th(t,state,model,x_eq,u_eq,ctrl)

% Select control
if ctrl == 1
    u = p_controller(t,state,x_eq,u_eq,model);
elseif ctrl == 2
    z = state(13:15);   % integral error states
    [u, z_dot] = pi_controller(t,state,x_eq,u_eq,model);
else 
    u = u_eq;
end


%-------------------%
% STATE DEFINITIONS %
%-------------------%
V_xb    = state(1);
V_yb    = state(2);
V_zb    = state(3); 
p       = state(4);
q       = state(5);
r       = state(6);
phi     = state(7);
theta   = state(8);
psi     = state(9);
x_f     = state(10);
y_f     = state(11);
z_f     = state(12);

% Gimbal check...
if abs(theta)*180/pi > 80
    t
    fprintf('Gimbal lock. Exiting...')
    return;
end

%-------------------%
% AERODYNAMIC MODEL %
%-------------------%
[g,W,hmat,Imat,F_b,M_b] = model(state,u);
Ixx = Imat(1,1);  Iyy = Imat(2,2);  Izz = Imat(3,3);
Ixy = -Imat(1,2); Ixz = -Imat(1,3); Iyz = -Imat(2,3);


%---------------------%
% EQUATIONS OF MOTION %
%---------------------%
pqr_term = [(Iyy-Izz)*q*r + Iyz*(q^2-r^2) + Ixz*p*q - Ixy*p*r; 
            (Izz-Ixx)*p*r + Ixz*(r^2-p^2) + Ixy*q*r - Iyz*p*q;
            (Ixx-Iyy)*p*q + Ixy*(p^2-q^2) + Iyz*p*r - Ixz*q*r];
rot_mat = [cos(theta)*cos(psi), sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi), cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
           cos(theta)*sin(psi), sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi), cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
           -sin(theta), sin(phi)*cos(theta), cos(phi)*cos(theta)];
euler_mat = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);
             0, cos(phi), -sin(phi);
             0, sin(phi)/cos(theta), cos(phi)/cos(theta)];

Vdot = g/W*F_b + g*[-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)] ...
     + [r*V_yb-q*V_zb; p*V_zb-r*V_xb; q*V_xb-p*V_yb];
pqrdot = Imat \ ( M_b + hmat*[p;q;r] + pqr_term );
xyzdot = rot_mat*[V_xb; V_yb; V_zb];
eulerdot = euler_mat*[p;q;r];
xdot_aircraft = [Vdot; pqrdot; eulerdot; xyzdot];

% Add integral states if necessary 
if ctrl == 2
    state_dot = [xdot_aircraft; z_dot];
else 
    state_dot = xdot_aircraft;
end 
