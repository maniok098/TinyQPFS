function sys_obs_dist = getObserver(omega0,D,Ts)
% disturbance observer
A = [0 1 0;
     0 0 1;
     0 -omega0^2 -2*D*omega0];
B = [0;
     0;
     omega0^2]; 
C = [1 0 0];

% % Parameter of Kalman filter
Q = 1e-3*diag([5,6,10]);
R = (5e-6/3)^2;
% solve Riccati equation to get covariance matrix P
[P,~,~] = icare(A',[],Q,[],[],[],-C'*R^(-1)*C);
% % Kalman gain
L = P*C'*R^(-1);

Ahat = A-L*C;
Bhat = [B, L];
Chat = eye(3);
Dhat = [0 0; 0 0; 0 0];
sys_obs_cont = ss(Ahat,Bhat,Chat,Dhat);
sys_obs_dist = c2d(sys_obs_cont,Ts);
end
