% SOLVE CLF - CBF OPTIMIZATION PROBLEM
% Uses CPLEX to solve
% minimize 0.5*u'*H(x)*u + k*delta^2
% s. t.    LfV(x) + LgV(x)u <= -gamma(V(x)) + delta
%          Lfh(x) + Lgh(x)u >= -alpha(h(x))
%
%
% By expressing it as:
% minimize 0.5*x'*H*x
% s. t.    Aineq*x <= bineq
%
% Sample use:
%   u = CLF_CBF_QP(x,Q,q,r);
%   [u, lambda] = CLF_CBF_QP(x,Q,q,r);
%
% INPUT: Elliptical obstacle parameters Q, q, r
%        Current state x
% OUTPUT: Optimal control input for above QP u
%         Correspdonding dual variables lambda
%[u(i,:) lambda] = Jankovic_CLF_CBF_QP(x(i-1,:),Q,q,r);

function [optimal_u, dual_variables, delta] = CLF_CBF_QP(x)
[Q, q, r] = get_ellipse_parameters();
sysdim = size(x,2);
k = 1.0;

H = [eye(sysdim) zeros(sysdim,1); zeros(1,sysdim) 2.0*k]; % construct H
Aineq = [LgV(x) -1; -Lgh(x) 0]; % construct Aineq
bineq = [-gain_gamma(V(x))-LfV(x); Lfh(x)+gain_alpha(h(x))]; % construct bineq
f = zeros(1,sysdim+1); % construct f (zero)
[u_opt, fval, exitflag, output, lambda] = cplexqp(H,f,Aineq,bineq);
dual_variables = lambda.ineqlin;
delta = u_opt(end);
optimal_u = u_opt(1:sysdim);
end

% Function defines PD matrix P for Lyapunov function
% V(x) = x'*P*x
function P = get_lyapunov_P()
%P = diag([1/4.5 1]);
%[P, q, r] = get_ellipse_parameters();
%P = P/r;
P = eye(2);
end

% Lyapunov function for system
% V(x) = x'*P*x
function lyapunov_V = V(x)
P = get_lyapunov_P();
lyapunov_V = x*P*x';
end


% VECTOR FUNCTIONS INPUT AND RETURN A ROW VECTOR
% Lie derivative of Lyapunov function along f
function Lie_fV = LfV(x)
[A, B] = get_linear_dynamics();
P = get_lyapunov_P();
Lie_fV = x*(P*A + A'*P)*(x');
end

% Lie derivative of Lyapunov function along g
function Lie_gV = LgV(x)
[A, B] = get_linear_dynamics();
P = get_lyapunov_P();
Lie_gV = 2*x*P*B;
end

% Safety function for forward invariance outside ellipse
function safety_h = h(x)
[Q, q, r] = get_ellipse_parameters();
safety_h = (x-q)*Q*(x-q)'-r^2;
end

% Lie derivative of safety function along f(x)
function Lie_fh = Lfh(x)
[A, B] = get_linear_dynamics();
[Q, q, r] = get_ellipse_parameters();
Lie_fh = 2*(x-q)*Q*A*(x');
end

% Lie derivative of safety function along g(x)
function Lie_gh = Lgh(x)
[Q, q, r] = get_ellipse_parameters();
[A, B] = get_linear_dynamics();
Lie_gh = 2*(x-q)*Q*B;
end

% Gain functions for optimization step

% Class K-infinity function alpha
% Characterizes gain on safety inequality
function a = gain_alpha(x)
a = x;
end

% Class K-infinity function gamma
% Characterizes gain on Lyapunov function inequality
function g = gain_gamma(x)
g = x;
end
