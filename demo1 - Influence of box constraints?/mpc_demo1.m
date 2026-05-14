clear
close all
clc
%% addpath
paths = {'../data','../utils'};
cellfun(@addpath, paths);

%% load trajectory and model
% get reference signal
load('traj_7pp.mat')
% time = traj_7pp.Time;
traj_KGT = traj_7pp.Variables;
traj_KGT = [traj_KGT;zeros(1000,4)];
Ts = 1/4e3; % sample time
time = Ts*(0:length(traj_KGT)-1)';

xdes = traj_KGT(:,1);
vdes = traj_KGT(:,2);
vd = vdes;
ades = traj_KGT(:,3);
jdes = traj_KGT(:,4);

% model parameters
% two-mass no fric
m1 = 140;
m2 = 470;
k = 6.0e7;
d = 1.0e4;
% friction
fc1 = 25;
fv1 = 300;
fc2 = 65;
fv2 = 380;
% transmission ratio
gamma = 0.04/(2*pi);
beta = 1e3; % smoothing factor of sign function
theta = [m1 m2 k d fc1 fv1 fc2 fv2 gamma beta];

% torque feedforward
lin_acc_2_torque = 3.902;  % nominal equivalent "mass"  % a_table to tau_m

clear m1 m2 k d fc1 fv1 fc2 fv2 gamma beta
%% Simulation with PPI and FF
x_sim = zeros(4,length(time));

% P position controller
Kv = 70;
% PI velocity controller
KpV = 155;
KiV = KpV*20*Ts;
% init u_PI
e_vm_previous = 0;
u_PI_previous = 0;

for k = 1:length(time)-1
  
  % read init condition
  xk = x_sim(:,k);
  
  % outer loop
  vm_des = vdes(k) + Kv*(xdes(k)-xk(2));
  
  % motor velocity error
  e_vm = vm_des - xk(3);
  % PI controller
  u_PI = u_PI_previous + KpV  *(e_vm - e_vm_previous) + KiV*e_vm;

  % control action
  uk = ades(k)*lin_acc_2_torque + u_PI;
  
  % one step ahead
  x_sim(:,k+1) = lpred_RK4(xk,uk,theta,Ts);
  
  % update for PI 
  e_vm_previous =  e_vm;
  u_PI_previous =  u_PI;

end

figure
subplot(211)
plot(time,xdes,time,x_sim(2,:))
ylabel('Position [m]')
title('Trajectory Tracking')
legend('x_{des}','x_{act}')

e_x = xdes'-x_sim(2,:);

subplot(212)
plot(time,e_x*1e6)
hold on
ylabel('Error [µm]')

%% setup model for MPC

% PT2I model
omega0 = 215.1584;
D0 = 0.44;
A = [0 1 0;
    0 0 1;
    0 -omega0^2 -2*D0*omega0];
B = [0 0 omega0^2]';
C = [1 0 0];

mdl_mpc = ss(A,B,C,0);
mdl_mpc_discrete = c2d(mdl_mpc,Ts);

Ad = mdl_mpc_discrete.A;
Bd = mdl_mpc_discrete.B;

% observer
sys_obs_dist = getObserver(omega0,D0,Ts);
Ad_obs = sys_obs_dist.A;
Bd_obs = sys_obs_dist.B;

clear A B C mdl_mpc mdl_mpc_discrete

%% MPC setup
% prediction horizon
N = 20;
% define weights
Q = diag([100 0 0])*1;   % stage cost of x(k)-xref(k)
R = 0.0001;             % stage cost of u(k)-uref(k)
P = 100*Q; % terminal cost of xN(N)-xref(N) % better choice
% set constraints
% constraints of x(k)
xlimits = [-0.8 0.8;   % position
        -1 1;       % velocity
        -20 20];    % acceleration
% constraints of u(k)
ulimits = [-2 2];   % commanded velocity
% terminal constraints of x(N) :
xlimitsEnd = xlimits;
% options for qpOASES
% options = qpOASES_options( 'maxIter',100 ); % -1 : default
options = qpOASES_options('MPC' ,'maxIter',100, 'printLevel', 0 ); % -1 : default
%% nominal MPC 

% init u_PI
e_vm_previous = 0;
u_PI_previous = 0;

% init sim data
x_sim_nMPC = zeros(4,length(time));
u_sim_nMPC = zeros(1,length(time));

% record computation time of QPs
runtime = 0;
vl_old = 0;

for k = 1:length(time)-N-1
  
  % read init condition
  xk = x_sim_nMPC(:,k);

  % inital state of system
  x0 = [xk(2);xk(4);(xk(4)-vl_old)/Ts];
  uref =   vd(k:k+N-1)';
  xref =   traj_KGT(k:k+N,1:3)';
  
  % set up optimization problem for quadprog
  if k ==1
      [H,f,~,~,Aeq,beq] = constructQP(Ad,Bd,x0,xref,uref,Q,R,P,xlimits,ulimits,xlimitsEnd);
      % setup lower and upper bounds for qpOASES
      lb = [repmat([xlimits(:,1);ulimits(:,1)],[N,1]);xlimitsEnd(:,1)];
      ub = [repmat([xlimits(:,2);ulimits(:,2)],[N,1]);xlimitsEnd(:,2)];
  else
    % update init condition and ref trajectory
     beq(1:3) = x0;
     f = [f(5:end-3); -Q*xref(:,end-1);-R*uref(:,end); -P*xref(:,end)];
  end
  
  tic
    if k ==1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
        [QP,solutionOL] = qpOASES_sequence('i', H,f,Aeq,lb,ub,beq,beq);
    else
        solutionOL = qpOASES_sequence('h',QP,f,lb,ub,beq,beq,options); % correct
    end
  runtime = runtime +toc;

  % read control action
  vm_des = solutionOL(4);
  
  % motor velocity error
  e_vm = vm_des - xk(3);
  % PI controller
  u_PI = u_PI_previous + KpV  *(e_vm - e_vm_previous) + KiV*e_vm;

  % apply control action
  x_next = lpred_RK4(xk,u_PI,theta,Ts);
  x_sim_nMPC(:,k+1) = x_next;
  u_sim_nMPC(:,k) = vm_des;
  
  % update for PI 
  e_vm_previous =  e_vm;
  u_PI_previous =  u_PI;
  
  % save the old table velocity for obtain acc for next step
  vl_old = xk(4);
end

e_x = xdes'-x_sim_nMPC(2,:);

figure(1)
subplot(212)
plot(time,e_x*1e6)

disp(['Calculation time ',num2str(runtime), 's'])


u_MPC_1 = u_sim_nMPC;
%% nominal MPC with only equality constraints

% init u_PI
e_vm_previous = 0;
u_PI_previous = 0;

% init sim data
x_sim_nMPC = zeros(4,length(time));
u_sim_nMPC = zeros(1,length(time));

% record computation time of QPs
runtime = 0;
vl_old = 0;

for k = 1:length(time)-N-1
  
  % read init condition
  xk = x_sim_nMPC(:,k);

  % inital state of system
  x0 = [xk(2);xk(4);(xk(4)-vl_old)/Ts];
  uref =   vd(k:k+N-1)';
  xref =   traj_KGT(k:k+N,1:3)';
  
  % set up optimization problem for quadprog
  if k ==1
      [H,f,~,~,Aeq,beq] = constructQP(Ad,Bd,x0,xref,uref,Q,R,P,xlimits,ulimits,xlimitsEnd);
      % Get rid of box constraints
      lb = [];
      ub = [];
  else
    % update init condition and ref trajectory
     beq(1:3) = x0;
     f = [f(5:end-3); -Q*xref(:,end-1);-R*uref(:,end); -P*xref(:,end)];
  end
  
  tic
    if k ==1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
        [QP,solutionOL] = qpOASES_sequence('i', H,f,Aeq,lb,ub,beq,beq);
    else
        solutionOL = qpOASES_sequence('h',QP,f,lb,ub,beq,beq,options); % correct
    end
  runtime = runtime +toc;

  % read control action
  vm_des = solutionOL(4);
  
  % motor velocity error
  e_vm = vm_des - xk(3);
  % PI controller
  u_PI = u_PI_previous + KpV  *(e_vm - e_vm_previous) + KiV*e_vm;

  % apply control action
  x_next = lpred_RK4(xk,u_PI,theta,Ts);
  x_sim_nMPC(:,k+1) = x_next;
  u_sim_nMPC(:,k) = vm_des;
  
  % update for PI 
  e_vm_previous =  e_vm;
  u_PI_previous =  u_PI;
  
  % save the old table velocity for obtain acc for next step
  vl_old = xk(4);
end
u_MPC_2 = u_sim_nMPC;

disp(['Calculation time ',num2str(runtime), 's'])

e_x = xdes'-x_sim_nMPC(2,:);

figure(1)
subplot(212)
plot(time,e_x*1e6)


figure(2)
subplot(211)
plot(time,u_MPC_1)
hold on
plot(time,u_MPC_2)
ylabel('u')

title('Box-Constrained MPC')
legend('BoxCon On','BoxCon Off',Location='best',NumColumns=2)

subplot(212)
plot(time,u_MPC_1-u_MPC_2)
xlabel('time')
ylabel('\Delta u')

%% Setup figs
figure(1)
subplot(211)
xlim([0,time(end)])
subplot(212)
xlim([0,time(end)])
legend('P','MPC1','MPC2',Location='best',NumColumns=3)

fig_size = [800 600];
fig_start = [200 -500];
figure(1)
set(gcf, 'Position', [fig_start fig_size])
figure(2)
set(gcf, 'Position', [fig_start fig_size])
%% rmpath
cellfun(@rmpath, paths);

