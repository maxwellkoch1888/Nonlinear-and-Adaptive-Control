function u = p_controller(state)
% omega_dot = f(x) + G(x)*controls

% Aerodynamic Properties
C_ell_beta = -0.0786;
C_ell_pbar = -0.3182;
C_ell_rbar = 0.0469;
C_ell_Lrbar = 0.1067;
C_ell_delta_a = -0.0741;
C_ell_delta_r = 0.0257;
C_m_0 = -0.0097;
C_m_alpha = 0.1766;
C_m_qbar = -4.8503;
C_m_delta_e = -0.5881;
C_n_beta = 0.2426;
C_n_pbar = 0.0131;
C_n_Lpbar = -0.1005;
C_n_rbar = -0.1787;
C_n_delta_a = -0.0276;
C_n_Ldelta_a = 0.0077;
C_n_delta_r = -0.0899;

S_w = 300;
b_w = 30;
cbar_w = 11.32;

C_L_0 = 0.0456;
C_L_alpha = 3.5791;

% Pull out actual rates
p = state(4); 
q = state(5); 
r = state(6); 
omega = state(4:6);

% Velocity and aerodynamic variables
V = sqrt(state(1)^2 + state(2)^2 + state(3)^2);
alpha = atan(state(3)/state(1));
beta = asin(state(2)/V);
pbar = p*b_w/(2*V);
qbar = q*cbar_w/(2*V);
rbar = r*b_w/(2*V);
C_L1 = C_L_0 + C_L_alpha*alpha;

% Density
[~,~,~,rho] = atmosisa(real(-state(12))*0.3048);  % meters
rho = rho * 0.001940320; % slug/ft^3

% Define desired angular acceleration with proportional controller
Kp = diag([2,2,2]);
omegadot_des = -Kp * (omega); % desired is zero 

% hmat
h_xb = 160;
h_yb = 0;
h_zb = 0;
hmat = [0, -h_zb, h_yb; h_zb, 0, -h_xb; -h_yb, h_xb, 0];

% Imat
Ixx = 9496;
Iyy = 55814;
Izz = 63100;
Ixy = 0;
Ixz = 982;
Iyz = 0;
Imat = [Ixx, -Ixy, -Ixz; -Ixy, Iyy, -Iyz; -Ixz, -Iyz, Izz];

pqr_term = [(Iyy-Izz)*q*r + Iyz*(q^2-r^2) + Ixz*p*q - Ixy*p*r; 
            (Izz-Ixx)*p*r + Ixz*(r^2-p^2) + Ixy*q*r - Iyz*p*q;
            (Ixx-Iyy)*p*q + Ixy*(p^2-q^2) + Iyz*p*r - Ixz*q*r];

% Dynamic inversion terms
Gamma = 1/2*rho*V^2*S_w * diag([b_w, cbar_w, b_w]);

C_states =  [C_ell_beta * beta + C_ell_pbar * pbar + (C_ell_Lrbar*C_L1 + C_ell_rbar) * rbar;
             C_m_0 + C_m_alpha*alpha + C_m_qbar*qbar;
             C_n_beta*beta + (C_n_Lpbar*C_L1 + C_n_pbar)*pbar + C_n_rbar*rbar];
C_control = [ C_ell_delta_a, 0, C_ell_delta_r;
             0, C_m_delta_e, 0;
             (C_n_Ldelta_a*C_L1 + C_n_delta_a), 0, C_n_delta_r ];

% Calculate f(x) 
f = Imat \ (Gamma * C_states + pqr_term + hmat*[p;q;r]);

% Calculate G(x)
G = Imat \ (Gamma * C_control);

% Dynamic inversion
delta = pinv(G)*(omegadot_des - f);

% Build full control vector
u(1:3) = delta;
u(4) = 0.282388376195832;

end
