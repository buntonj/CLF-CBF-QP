[Q, q, r] = get_ellipse_parameters();
% PLOT ELLIPTICAL OBSTACLE
L = chol(Q);
ellipse_pts = linspace(0,2*pi,300);
z = [cos(ellipse_pts); sin(ellipse_pts)]*r; % creates circle of correct radius
R = rotation_matrix(pi/4);
ellipse = L^-1 * z + q';    % stretch and translate circle appropriately
clf;

% PLOT GENERATED TRAJECTORY and VECTOR FIELD
limits = [min(xrange) max(xrange) min(yrange) max(yrange)];
axis(limits)
hold on
plot(ellipse(1,:),ellipse(2,:),'black','LineWidth',3)
quiver(plotx,ploty,plotux,plotuy,'red')
for i = 1:numtrajectories
    scatter(x(1,1,i),x(1,2,i),50,'blue','filled')
    plot(x(:,1,i),x(:,2,i),'blue','LineWidth',2)
end
legend('Obstacle Boundary','Input Vector','Initial State','Trajectory')