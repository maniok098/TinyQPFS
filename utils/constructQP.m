function [H,f,Aineq,bineq,Aeq,beq] = constructQP(Ad,Bd,x0,xref,uref,Q,R,P,xlim,ulim,xlimEnd)
% formulate linear MPC in quadratic programming 
%  - Adapted form for Matlab Quadprog 
%  - Author: Haijia Xu, 20.02.2023

% Prerequisite: PT2I approximation of PI-controlled inner velocity loop
% .                | 0  1  0 |     | 0  |
% x = A*x + B*u  = | 0  0  1 | x + | 0  | u
%                  | a0 a1 a2|     | b1 |
% 

% Linear system dynamics (discrete) :
%  x(k+1) = Ad*x(k) + Bd*u(k)

% Optimal control problem of linear MPC:
% Cost function : 
%                 J = 0.5*|x(N)-xref(N)|^2_P 
%                   + 0.5*sum_k^{N-1}  |x(k)-xref(k)|^2_Q +  |u(k)-uref(k)|^2_R 
% subject to    :
%                 x(k+1) = Ad*x(k) + Bd*u(k),  x(0) = x0
%                 xlim : nStates x 2  [min, max] for x(k), k=0:N-1
%                 ulim : nInputs x 2  [min, max] for u(k), k=0:N-1
%              xlimEnd : nStates x 2  [min, max] for x(N) (terminal constraint)
% 
% Prediction horizon is determined by dimension of xref and uref
%  xref: reference of [x0, x1, ..., xN-1, xN]   dimension : nStates x N+1
%  uref: reference of [u0, u1, ..., uN-1] dimension : nInputs x N
%     Q: stage cost of x
%     R: stage cost of u
%     P: terminal cost of x

% Target QP: min  0.5x^THx +f^T*x
%            s.t. Ax    <= b
%                 Aeq x  = beq
%  for Quadprog

% read dimensions
nStates = length(Q);
nInputs = length(R);
N = size(uref,2); % prediction horizon

%% Cost function
% Construct H 
% Target : diagonal of Q R Q R... (N times) then P
%  |Q       |
%  | R      |
%  |  Q     |
%  |   R    |
%  |    ... |
%  |       P|
% Dimension: square matrix -> N*(nStates+nInputs)+nInputs
length_opt_var = N*(nInputs+nStates)+nStates;

% repeat diagonal matrix using Kronecker tensor product
H = blkdiag(kron(eye(N),blkdiag(Q,R)),P); 

% Construct f
%  |-Q*xref(0)  |
%  |-R*uref(0)  |
%  |...         |
%  |-Q*xref(N-1)|
%  |-R*uref(N-1)|
%  |-P*xref(N)  |
f = zeros(length_opt_var,1);
f(1:length_opt_var-nStates) = reshape([-Q*xref(:,1:N);-R*uref],[length_opt_var-nStates,1]);
f(length_opt_var-nStates+1:end) = -P*xref(:,N+1);

%% Equality constraints
% Construct Aeq
%  |I                |
%  |A B -I           |
%  |     A B -I      |
%  |          A B -I |
%  |              ...|
% Dimension : 
%  - row:    (N+1)*nStates
%  - column: N*(nStates+nInputs)+nStates
Aeq =[ [eye(nStates) , zeros(nStates,N*(nInputs+nStates))] ;
       [kron(eye(N),[Ad,Bd]),zeros(N*nStates,nStates)] + ...
       [zeros(N*nStates,nStates),kron(eye(N),[0*Bd,-eye(nStates)])]  ];

% Construct beq
%  | x0|
%  | 0 |
%  | 0 |
%  |...|
beq = [x0(:);zeros(N*nStates,1)];

%% Inequality constraints
% Construct Aineq
%  |-I_x         |
%  | I_x         |
%  |    -I_u     |
%  |     I_u     |
%  |     ...     |
%  |        -I_x |
%  |         I_x |
% Dimension : 
%  - row :    2*N*(nStates+nInputs) + 2*nStates
%  - column : N*(nStates+nInputs)+nStates
tmp1  = [-eye(nStates);eye(nStates)];
tmp2  = [-eye(nInputs);eye(nInputs)];
Aineq = blkdiag(kron(eye(N),blkdiag(tmp1,tmp2)),tmp1);

% Construct bineq
%  |-xmin(0)   |
%  | xmax(0)   |
%  |-umin(0)   |
%  | umax(0)   |
%  |     ...   |
%  |-xmin(N-1) |
%  | xmax(N-1) |
%  |-umin(N-1) |
%  | umax(N-1) |
%  |-xmin(N)   |
%  | xmax(N)   |
%  Notice: - Constants for 0:N-1 are assumed to be the same
%          - Terminal constaints of x(N) can be different
bineq = [repmat([-xlim(:,1);xlim(:,2);-ulim(:,1);ulim(:,2)],[N,1]); -xlimEnd(:,1);xlimEnd(:,2)];

end
