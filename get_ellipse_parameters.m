function [Q, q, r] = get_ellipse_parameters()
Q = diag([1/9 1]);
R = rotation_matrix(0);
%Q = R'*[1 -0.5; -0.5 1]*R;
q = [0 2];
r = 1;
end