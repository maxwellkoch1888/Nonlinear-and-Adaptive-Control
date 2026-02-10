clear all; close all; 

syms x1 x2 

f1 = -4*x1 + 2*x2 - 3 * x1^2 - x1^3; 
f2 = 2*x1 - 4*x2 - 3*x2^2 - x2^3; 

A = [diff(f1,x1), diff(f1,x2); diff(f2,x1), diff(f2,x2)]; 
A = double(subs(A, {x1,x2}, {0,0}));


% Lyapunov 
P = lyap(A.', eye(2)); % dot does element wise 
V = simplify([x1,x2] * P * [x1;x2]);
Vd = simplify([f1,f2] * P * [x1;x2] + [x1,x2]* P* [f1;f2]); % derivative of V

% Numerical estimates
c = 0.48; % try different c values to make Vd_max less than or equal to zero, ideally less than
x1 = linspace(-3,3); 
x2 = x1; 
[x1,x2] = meshgrid(x1,x2); 
V_c = eval(vectorize(V)); 
Vd_c = eval(vectorize(Vd)); 

% Is the set V <= c contained in Vd <= 0?  Vd_max shoudl be non-positive. 
id = V_c <= c; 
Vd_max = max(Vd_c(id))

% Make sure we are still in the domain
RoA = zeros(size(x1)); 
RoA(id) = 1; 
figure, contourf(x1,x2,RoA), colorbar 