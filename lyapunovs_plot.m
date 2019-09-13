num_lyapunovs = 100;
num_points = 300;

lyapunov = zeros(num_lyapunovs, num_points);
ellipse_pts = linspace(0,2*pi,num_points);
[Q,q,r] = get_ellipse_parameters();
L = chol(Q);
z = [cos(ellipse_pts); sin(ellipse_pts)]*r; % creates circle of correct radius
ellipse = L^-1 * z + q';
ellipse = [ellipse; zeros(1,300)];

P = zeros(num_lyapunovs,2,2);
P(:,1,1) = linspace(1/9,1,num_lyapunovs);
P(:,2,2) = 1;
[A,B] = get_linear_dynamics();
problem = problem_instance(A,B,squeeze(P(1,:,:)),Q,q,r);
problem.k = 1.0;
max_points = zeros(num_lyapunovs,2,3);
deltas = zeros(num_lyapunovs,num_points);
for j = 1:num_lyapunovs
    problem.P = squeeze(P(j,:,:));
    for i = 1:num_points
        lyapunovs(j,i) = problem.V(ellipse(1:2,i)');
        [u, l, deltas(j,i)] = problem.CLF_CBF_QP(ellipse(1:2,i)');
    end
    
    max_delta(j) = max(deltas(j,:),[],2);
    
    max_val1(j) = max(lyapunovs(j,:),[],2);
    max_idx1 = find(lyapunovs(j,:)==max_val1(j));
    num_max = size(max_idx1,2); 
    max_val2(j) = max(lyapunovs(j,(lyapunovs(j,:) < max_val1(j))));
    max_idx2 = find(lyapunovs(j,:) == max_val2(j));
    
    max_points(j,:,:) = [ellipse(1,max_idx1), ellipse(2,max_idx1), max_val1(j);
                  ellipse(1,max_idx2), ellipse(2,max_idx2), max_val2(j)];
end

%ellipse = plot3(ellipse(1,:),ellipse(2,:),ellipse(3,:),'LineWidth',2);
%axis(limits);
%grid on;
%hold on;
%lyap_plot = plot3(ellipse(1,mm:),ellipse(2,:),lyapunovs(1,:),'LineWidth',2);
 
limits = [-4 4 0 4 0 max(lyapunovs,[],'all')];
h = figure;
axis(limits);
filename = 'lyapunovs.gif';
for j = num_lyapunovs:-1:1
    ellipse_plot = plot3(ellipse(1,:),ellipse(2,:),ellipse(3,:),'black','LineWidth',2);
    axis(limits);
    grid on;
    hold on;
    %delete(lyap_plot)
    plot3(ellipse(1,:),ellipse(2,:),lyapunovs(j,:),'blue','LineWidth',2);
    scatter3(max_points(j,:,1),max_points(j,:,2),max_points(j,:,3),75,'filled','red');
    title_string = sprintf("P = diag( %f, %f)", round(P(j,1,1),3),round(P(j,2,2),3));
    title(title_string)
    zlabel('V(x)');
    
    drawnow;
    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
    
    if j == num_lyapunovs
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append','DelayTime',0.02);
    end
    %pause(0.01);
    hold off;
end
for j = 1:num_lyapunovs
    ellipse_plot = plot3(ellipse(1,:),ellipse(2,:),ellipse(3,:),'black','LineWidth',2);
    axis(limits);
    grid on;
    hold on;
    %delete(lyap_plot)
    plot3(ellipse(1,:),ellipse(2,:),lyapunovs(j,:),'blue','LineWidth',2);
    scatter3(max_points(j,:,1),max_points(j,:,2),max_points(j,:,3),75,'filled','red');
    title_string = sprintf("P = diag( %f, %f)", round(P(j,1,1),3),round(P(j,2,2),3));
    title(title_string)
    zlabel('V(x)');
    
    drawnow;
    frame = getframe(h);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im,256);
    %pause(0.01);
    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append','DelayTime',0.02);
    
    hold off;
end