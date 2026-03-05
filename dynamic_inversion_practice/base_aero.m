function [g,W,hmat,Imat,F_b,M_b] = base_aero(state,control)

% State Variables
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

% Derived Quantities
V       = sqrt(V_xb^2+V_yb^2+V_zb^2);
alpha   = atan(V_zb/V_xb);
beta    = asin(V_yb/V);

% Control Variables
delta_a = control(1);
delta_e = control(2);
delta_r = control(3);
tau     = control(4);

% Thrust Model (See Section II.D)
H = -z_f;
[~,~,~,rho] = atmosisa(real(H)*0.3048);   rho = rho*0.001940320;
[~,~,~,rho_0] = atmosisa(0);        rho_0 = rho_0*0.001940320;

% For a:
a_idle  = 1.0104 + 2.9484e-5*H - 3.8270e-10*H^2;
a_mil   = 1.0148 + 3.1355e-5*H - 4.2106e-10*H^2;
a_max   = 1.0225 + 3.1984e-5*H - 4.3617e-10*H^2;

% For T0:
T0_idle =  3145 - 0.4185*H + 1.8313e-5*H^2;
T0_mil  = 11716 + 0.1156*H + 0.3474e-5*H^2;
T0_max  = 20341 + 0.1454*H + 0.9283e-5*H^2;

% For T1:
T1_idle = -4.3491 - 4.9703e-4*H + 1.3557e-8*H^2;
T1_mil  =  3.5689 + 0.1409e-4*H - 0.3982e-8*H^2;
T1_max  =  1.9886 + 6.3926e-4*H - 2.4428e-8*H^2;

% For T2:
T2_idle = -0.2321e-3 + 5.5629e-7*H - 2.0550e-11*H^2;
T2_mil  = -3.9793e-3 + 2.6931e-7*H + 0.5281e-11*H^2;
T2_max  =  3.5201e-3 + 0.7574e-7*H + 2.6665e-11*H^2;

% For T_{}:
T_idle = (rho/rho_0)^a_idle * ( T0_idle + T1_idle*V + T2_idle*V^2 );
T_mil  = (rho/rho_0)^a_mil  * ( T0_mil  + T1_mil *V + T2_mil *V^2 );
T_max  = (rho/rho_0)^a_max  * ( T0_max  + T1_max *V + T2_max *V^2 );

% Fot T:
if tau <= 0.77
    P1 = 64.94*tau;
else
    P1 = 217.38*tau - 117.38;
end

if P1 < 50
    T = T_idle + (T_mil-T_idle)*P1/50;
else
    T = T_mil + (T_max-T_mil)*(P1-50)/50;
end

% Geometric properties (See Table A.1)
S_w = 300;
b_w = 30;
cbar_w = 11.32;
h_xb = 160;
h_yb = 0;
h_zb = 0;
hmat = [0, -h_zb, h_yb; h_zb, 0, -h_xb; -h_yb, h_xb, 0];

% Weight and inertia (See Table A.2)
g = 32.1278;
W = 20500;
Ixx = 9496;
Iyy = 55814;
Izz = 63100;
Ixy = 0;
Ixz = 982;
Iyz = 0;
Imat = [Ixx, -Ixy, -Ixz; -Ixy, Iyy, -Iyz; -Ixz, -Iyz, Izz];

% Aero Force Coefficients (See Table A.4)
C_L_0 = 0.0456;
C_L_alpha = 3.5791;
C_L_qbar = 3.3916;
C_L_delta_e = 0.5652;
C_S_beta = -0.9009;
C_S_pbar = -0.0153;
C_S_Lpbar = 0.3318;
C_S_rbar = 0.4357;
C_S_delta_a = 0.0656;
C_S_delta_r = 0.1698;
C_D_0 = 0.0218;
C_D_L = -0.0340;
C_D_L2 = 0.1834;
C_D_S2 = 0.6081;
C_D_Spbar = 0.0768;
C_D_qbar = 0.0368;
C_D_Lqbar = 0.7750;
C_D_L2qbar = -0.1844;
C_D_Srbar = -0.7239;
C_D_delta_e = -0.0032;
C_D_Ldelta_e = 0.1775;
C_D_delta_e2 = 0.2854;
C_D_Sdelta_a = 0.1118;
C_D_Sdelta_r = 0.1604;

% Aero Force Coefficients (See Eqs. 8-10)
C_L_1 = C_L_0 + C_L_alpha*alpha;
C_S_1 = C_S_beta*beta;
pbar = p*b_w/(2*V);
qbar = q*cbar_w/(2*V);
rbar = r*b_w/(2*V);
C_L = C_L_1 + C_L_qbar*qbar + C_L_delta_e*delta_e;
C_S = C_S_beta*beta + (C_S_Lpbar*C_L_1+C_S_pbar)*pbar + ...
      C_S_rbar*rbar + C_S_delta_a*delta_a + C_S_delta_r*delta_r;
C_D = C_D_0 + C_D_L*C_L_1 + C_D_L2*C_L_1^2 + C_D_S2*C_S_1^2 + ...
      C_D_Spbar*C_S_1*pbar + (C_D_L2qbar*C_L_1^2 + C_D_Lqbar*C_L_1 + C_D_qbar)*qbar + ...
      C_D_Srbar*C_S_1*rbar + C_D_Sdelta_a*C_S_1*delta_a + ...
      (C_D_Ldelta_e*C_L_1 + C_D_delta_e)*delta_e + C_D_delta_e2*delta_e^2 + ...
      C_D_Sdelta_r*C_S_1*delta_r;
  
  
% Aero Moment Coefficients (See Table A.5)
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

% Aero Moment Coefficients (See Eqs. 11-13)
C_ell = C_ell_beta*beta + C_ell_pbar*pbar + ...
        (C_ell_Lrbar*C_L_1 + C_ell_rbar)*rbar + C_ell_delta_a*delta_a + C_ell_delta_r*delta_r;
C_m = C_m_0 + C_m_alpha*alpha + C_m_qbar*qbar + C_m_delta_e*delta_e;
C_n = C_n_beta*beta + (C_n_Lpbar*C_L_1+C_n_pbar)*pbar + C_n_rbar*rbar + ...
      (C_n_Ldelta_a*C_L_1 + C_n_delta_a)*delta_a + C_n_delta_r*delta_r;

% Resultant Body Forces (See Eq. 38)
F_b = 1/2*rho*V^2*S_w * [C_L*sin(alpha) - C_S*cos(alpha)*sin(beta) - C_D*cos(alpha)*cos(beta);
     C_S*cos(beta) - C_D*sin(beta);
     -C_L*cos(alpha) - C_S*sin(alpha)*sin(beta) - C_D*sin(alpha)*cos(beta)] + ...
     [T;0;0];
    
% Resultant Body Moments (See Eq. 39)
M_b = 1/2*rho*V^2*S_w*[b_w*C_ell; cbar_w*C_m; b_w*C_n];