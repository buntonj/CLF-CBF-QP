% JANKOVIC's MODIFIED CLF - CBF OPTIMIZATION PROBLEM
% Uses CPLEX to solve
% minimize 0.5*u'*H*u + k*delta'*delta
% s. t.    LfV(x) + LgV(x)u <= -gamma(V(x)) - LgV(x)*delta
%          Lfh(x) + Lgh(x)u >= -alpha(h(x))    
%
% Sample use:
%   u = Jankovic_CLF_CBF_QP(x,Q,q,r);
%   [u, lambda] = Jankovic_CLF_CBF_QP(x,Q,q,r);
%
% INPUT: Elliptical obstacle parameters Q, q, r
%        Current state x
% OUTPUT: Optimal control input for above QP u
%         Correspdonding dual variables lambda
%[u(i,:) lambda] = Jankovic_CLF_CBF_QP(x(i-1,:),Q,q,r);

function [optimal_u, dualvars, delta] = Jankovic_CLF_CBF_QP(x,Q,q,r)
sysdim = size(x,2);
k = 1.0;

H = blkdiag(eye(sysdim), eye(sysdim)*k); % construct H

Aineq = [LgV(x) LgV(x); 
        -Lgh(x,Q,q,r) zeros(1,sysdim)]; % construct Aineq

bineq = [-jankovic_gamma(V(x)+LfV(x)); 
        Lfh(x,Q,q,r)+gain_alpha(h(x,Q,q,r))]; % construct bineq

f = zeros(1,sysdim*2); % construct f (zero)

[u_opt, fval, exitflag, output, lambda] = cplexqp(H,f,Aineq,bineq);
dualvars = lambda.ineqlin;
delta = u_opt(sysdim+1:end);
optimal_u = u_opt(1:sysdim);
end

% Function defines PD matrix P for Lyapunov function
% V(x) = x'*P*x
function P = get_lyapunov_P()
P = diag([1 1]);
end

% Lyapunov function for system
% V(x) = x'*x
function lyapunov_V = V(x)
P = get_lyapunov_P();
lyapunov_V = x*P*x';
end


% VECTOR FUNCTIONS INPUT AND RETURN A ROW VECTOR
% Lie derivative of Lyapunov function along f
function Lie_fV = LfV(x)
Lie_fV = 0;
end

% Lie derivative of Lyapunov function along g
function Lie_gV = LgV(x)
P = get_lyapunov_P();
Lie_gV = 2*x*P;
end

% Safety function for forward invariance outside ellipse
function safety_h = h(x,Q,q,r)
safety_h = (x-q)*Q*(x-q)'-r^2;
end

function grad_h = h_gradient(x,Q,q,r)

end

% Lie derivative of safety function along f(x)
function Lie_fh = Lfh(x,Q,q,r)
Lie_fh = 0;
end

% Lie derivative of safety function along g(x)
function Lie_gh = Lgh(x,Q,q,r)
Lie_gh = 2*(x-q)*Q;
end

% Class K-infinity function alpha
% Characterizes gain on safety inequality
function a = gain_alpha(x)
a = x;
end

% Gain function gamma defined in Jankovic (2018)
% applies a positive, >= 1 gain if argument is positive
% applies no gain otherwise
function gamma = jankovic_gamma(x)
g = 1.0; % gain on positive inputs
if (x < 0)
    gamma = x;
else
    gamma = g*x;
end
end
