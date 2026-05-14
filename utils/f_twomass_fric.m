
function x_dot = f_twomass_fric(x,u,theta)
% dynamics of two-mass model with friction
% xxu 23.02.2023

% x : [x_m x_T v_m v_T]^T
% u : tau_m

% read parameters 
m1 = theta(1);
m2 = theta(2);
k = theta(3);
d = theta(4);
% friction
fc1 = theta(5);
fv1 = theta(6);
fc2 = theta(7);
fv2 = theta(8);

% fixed transmission ratio
% gamma = 0.04/(2*pi); 
gamma = theta(9);

beta = theta(10);

% read states and input
xm = x(1);
xT = x(2);
vm = x(3);
vT = x(4);

tau_m = u(1);

% compute forces

% friction
% beta = 1e3; % smooting factor to replace sign(v)
fric_m = fc1*tanh(beta*vm) + fv1* vm;
fric_T = fc2*tanh(beta*vT) + fv2* vT;

% force of spring-damper
Fn = k*(xT - xm) + d*(vT - vm);

% dynamics
dx1 = vm;
dx2 = vT;
dx3 = 1/m1* (tau_m/gamma + Fn - fric_m );
dx4 = 1/m2* (-Fn - fric_T);

x_dot = [dx1; 
         dx2;
         dx3;
         dx4];

end