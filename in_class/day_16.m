clear all; close all; 

%% Choose beta = 1

% Generate numerical results
f = @(t,x) [x(2); -4*x(1) - 2*x(2) + x(2)^3]; 

for i = 1:1000; 
    x0 = randn(2,1); 
    X0(:,i) = x0; 
    [t,x] = ode45(f,[0,10], x0); 
    xF(i) = norm(x(end,:)); 
end 

id = xF < 0.1; 

figure, plot(X0(1,id), X0(2,id), 'g+'), grid on, hold on 
plot(X0(1,~id),X0(2,~id),'r+')
xlabel('x1'), ylabel('x2')

theta = linspace(0,2*pi);
rx = cos(theta);
ry = sin(theta); 

% Numerical Estimate 1
Q = eye(2); 
c = 0.3; % if Vd_max > 0, decrease c
[Vd_max,P] = RoA_Estimator(Q,c)
R = inv(sqrtm(P/c))*[rx;ry];
hold on, plot(R(1,:),R(2,:),'r-', 'linewidth',2)

% Numerical Estimate 2
Q = [1,0;0,100]; 
c = 51; % if Vd_max > 0, decrease c
[Vd_max,P] = RoA_Estimator(Q,c)
R = inv(sqrtm(P/c))*[rx;ry];
hold on, plot(R(1,:),R(2,:),'r--', 'linewidth',2)

% Numerical Estimate 3
Q = 2*randn(2,2); Q = Q'*Q; 
c = 0.009; % if Vd_max > 0, decrease c
[Vd_max,P] = RoA_Estimator(Q,c)
R = inv(sqrtm(P/c))*[rx;ry];
hold on, plot(R(1,:),R(2,:),'r-', 'linewidth',2)

function[Vd_max,P] = RoA_Estimator(Q,c)
    syms x1 x2 
    
    f1 = x2; 
    f2 = -4*x1-2*x2+x2^3; 
    
    A = [diff(f1,x1), diff(f1,x2); diff(f2,x1), diff(f2,x2)]; 
    A = double(subs(A, {x1,x2}, {0,0}));
    
    
    % Lyapunov 
    P = lyap(A.', Q); % dot does element wise 
    V = simplify([x1,x2] * P * [x1;x2]);
    Vd = simplify([f1,f2] * P * [x1;x2] + [x1,x2]* P* [f1;f2]); % derivative of V
    
    % Numerical estimates
    x1 = linspace(-1,1); 
    x2 = linspace(-2,2); 
    [x1,x2] = meshgrid(x1,x2); 
    V_c = eval(vectorize(V)); 
    Vd_c = eval(vectorize(Vd)); 
    
    % Is the set V <= c contained in Vd <= 0?  Vd_max shoudl be non-positive. 
    id = V_c <= c; 
    Vd_max = max(Vd_c(id));
end 
