[Q,q,r] = get_ellipse_parameters();
[A,B] = get_linear_dynamics();
p = problem_instance(A,B,eye(2),Q,q,r);
x = linspace(-4,4,500);
y = linspace(0,4,500);
[x,y] = meshgrid(x,y);
V = x;
for i = 1:500
for j = 1:500
state = [x(i,j) y(i,j)];
if not(p.in_obstacle(state))
V(i,j) = p.V(state);
else
V(i,j) = NaN;
end
end
end
surf(x,y,V,'edgecolor','none')