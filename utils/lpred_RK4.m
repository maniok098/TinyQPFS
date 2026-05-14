function x_next = lpred_RK4(xk,uk,theta,Ts)
% xk: x(k)
% uk: u(k)
% theta: model parameters
% Ts: sampling time

% selection dynamics function
f_fun = @f_twomass_fric;
nSteps = 10;
Ts = Ts/nSteps;
    for i = (1:nSteps)
        s1 = f_fun(xk,uk,theta);
        s2 = f_fun(xk+0.5*Ts*s1,uk,theta);
        s3 = f_fun(xk+0.5*Ts*s2,uk,theta);
        s4 = f_fun(xk+s3*Ts,uk,theta);
        
        % predicted x(k+1)
        x_next = xk + (1/6)*(s1+2*s2+2*s3+s4)*Ts;
        xk = x_next;
    end
end