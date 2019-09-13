classdef problem_instance < handle
    properties
        A
        B
        k = 1.0
        jankovic_gain = 1.0
        P
        Q
        q
        r
        sysdim
        ax
        ay
        ux
        uy
        limits
        x
        numtrajectories = 0
    end
    methods
        function self = problem_instance(A,B,P,Q,q,r)
            self.A = A;
            self.B = B;
            self.P = P;
            self.Q = Q;
            self.q = q;
            self.r = r;
            self.sysdim = size(A,2);
        end
        
        % SOLVE CLF - CBF OPTIMIZATION PROBLEM
        % Uses CPLEX to solve
        % minimize 0.5*u'*H(x)*u + k*delta^2
        % s. t.    LfV(x) + LgV(x)u <= -gamma(V(x)) + delta
        %          Lfh(x) + Lgh(x)u >= -alpha(h(x))
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
        function [optimal_u, dual_variables, delta] = CLF_CBF_QP(self,x)
            H = [eye(self.sysdim) zeros(self.sysdim,1); zeros(1,self.sysdim) 2.0*self.k]; % construct H
            Aineq = [self.LgV(x) -1; -self.Lgh(x) 0]; % construct Aineq
            bineq = [-self.gain_gamma(self.V(x))-self.LfV(x); self.Lfh(x)+self.gain_alpha(self.h(x))]; % construct bineq
            f = zeros(1,self.sysdim+1); % construct f (zero)
            [u_opt, ~, ~, ~, lambda] = cplexqp(H,f,Aineq,bineq);
            dual_variables = lambda.ineqlin;
            delta = u_opt(end);
            optimal_u = u_opt(1:self.sysdim);
        end
        
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
        function [optimal_u, dualvars, delta] = Jankovic_CLF_CBF_QP(self,x)
        H = blkdiag(eye(self.sysdim), eye(self.sysdim)*self.k); % construct H
        Aineq = [self.LgV(x) self.LgV(x); 
                -self.Lgh(x) zeros(1,self.sysdim)]; % construct Aineq
        bineq = [-self.jankovic_gamma(self.V(x)+self.LfV(x)); 
                self.Lfh(x)+self.gain_alpha(self.h(x))]; % construct bineq
        f = zeros(1,self.sysdim*2); % construct f (zero)
        [u_opt, ~, ~, ~, lambda] = cplexqp(H,f,Aineq,bineq);
        dualvars = lambda.ineqlin;
        delta = u_opt(self.sysdim+1:end);
        optimal_u = u_opt(1:self.sysdim);
        end

        % Differential equations governing motion
        function odefcn = xdot(self,t,x,u)
            odefcn = self.A*x + self.B*(u);
        end
        
        % Lyapunov function for system
        % V(x) = x'*P*x
        function lyapunov_V = V(self,x)
            lyapunov_V = x*self.P*x';
        end

        % VECTOR FUNCTIONS INPUT AND RETURN A ROW VECTOR
        % Lie derivative of Lyapunov function along f
        function Lie_fV = LfV(self,x)
            Lie_fV = x*(self.P*self.A + (self.A)'*self.P)*(x');
        end

        % Lie derivative of Lyapunov function along g
        function Lie_gV = LgV(self,x)
            Lie_gV = 2*x*self.P*self.B;
        end

        % Safety function for forward invariance outside ellipse
        function safety_h = h(self,x)
            safety_h = (x-self.q)*self.Q*(x-self.q)'-self.r^2;
        end

        function obst_test = in_obstacle(self, x)
            %state = x - self.q;
            obst_test = ((x-self.q)*self.Q*(x-self.q)' < self.r^2);
        end
        
        % Lie derivative of safety function along f(x)
        function Lie_fh = Lfh(self,x)
            Lie_fh = 2*(x-self.q)*self.Q*self.A*(x');
        end

        % Lie derivative of safety function along g(x)
        function Lie_gh = Lgh(self,x)
            Lie_gh = 2*(x-self.q)*self.Q*self.B;
        end

        % Class K-infinity function alpha
        % Characterizes gain on safety inequality
        function a = gain_alpha(self,x)
            a = x;
        end

        % Class K-infinity function gamma
        % Characterizes gain on Lyapunov function inequality
        function g = gain_gamma(self,x)
            g = x;
        end
        
        % Gain function gamma defined in Jankovic (2018)
        % applies a positive, >= 1 gain if argument is positive
        % applies no gain otherwise
        function gamma = jankovic_gamma(self,x)
            if (x < 0)
                gamma = x;
            else
                gamma = self.jankovic_gain*x;
            end
        end
        
        function generate_trajectories(self,method,x0,tsteps)
        % GENERATE TRAJECTORY
        %x0 = [x0_1 ; x0_2 ;...]
        % each row of x0 is a starting point for x
        self.numtrajectories = size(x0,1);
        stepsize = 1E-3; % step size between iterations
        %u = zeros(tsteps,sysdim,numtrajectories); % create vector of inputs
        self.x = zeros(tsteps,self.sysdim,self.numtrajectories); % create vector of states
        t = linspace(0,stepsize*tsteps,tsteps); % generate vector of times
        self.x(1,:,:) = x0;
        
        if method == 'standard'
            for j = 1:self.numtrajectories
                for i = 2:tsteps
                    % SOLVE STANDARD CLF - CBF OPTIMIZATION PROBLEM
                    u = self.CLF_CBF_QP(self.x(i-1,:,j));

                    % PERFORM INTEGRATION STEP
                    % Using ode45 for laziness
                    tspan = [t(i-1),t(i)]; % time interval for current iteration
                    [~, sol_x] = ode45(@(t,y) self.xdot(t,y,u), tspan, self.x(i-1,:,j)); % use ode45 to integrate
                    self.x(i,:,j) = sol_x(end,:);
                end
            end
        elseif method == 'jankovic'
           for j = 1:self.numtrajectories
                for i = 2:tsteps
                    % SOLVE STANDARD CLF - CBF OPTIMIZATION PROBLEM
                    u = self.Jankovic_CLF_CBF_QP(self.x(i-1,:,j));

                    % PERFORM INTEGRATION STEP
                    % Using ode45 for laziness
                    tspan = [t(i-1),t(i)]; % time interval for current iteration
                    [~, sol_x] = ode45(@(t,y) self.xdot(t,y,u), tspan, self.x(i-1,:,j)); % use ode45 to integrate
                    self.x(i,:,j) = sol_x(end,:);
                end
            end 
        end
        end
        
        function vector_field(self, method, input_limits)
            % GENERATE VECTOR FIELD FOR INPUTS
            % limits = [xmin xmax ymin ymax]
            % method is a STRING 'standard' or 'jankovic'
            self.limits = input_limits;
            xrange = linspace(self.limits(1),self.limits(2),100);
            yrange = [self.limits(3):(xrange(2)-xrange(1)):self.limits(4)];
            %xrange = [-4:0.1:4];
            %yrange = [0:0.1:8];
            [x, y] = meshgrid(xrange,yrange);
            self.ux = x;
            self.uy = y;
            self.ax = x;
            self.ay = y;

            if method == 'standard'
                for i = 1:size(self.ax,1)
                    for j = 1:size(self.ax,2)
                        state = [self.ax(i,j) self.ay(i,j)];
                        if state == self.q
                            self.ux(i,j) = 0;
                            self.uy(i,j) = 0;
                            self.ax(i,j) = 0;
                            self.ay(i,j) = 0;
                            continue
                        end

                        u_opt = self.CLF_CBF_QP(state);
                        self.ux(i,j) = u_opt(1)/sqrt(u_opt(1)^2 + u_opt(2)^2);
                        self.uy(i,j) = u_opt(2)/sqrt(u_opt(1)^2 + u_opt(2)^2);
                        %for not-normalized vector field
                        %plotux(i,j) = u_opt(1);
                        %plotuy(i,j) = u_opt(2);
                    end
                end
            elseif method == 'jankovic'
                for i = 1:size(self.ax,1)
                    for j = 1:size(self.ax,2)
                        state = [self.ax(i,j) self.ay(i,j)];
                        if state == self.q
                            self.ux(i,j) = 0;
                            self.uy(i,j) = 0;
                            self.ax(i,j) = 0;
                            self.ay(i,j) = 0;
                            continue
                        end

                        u_opt = self.Jankovic_CLF_CBF_QP(state);
                        self.ux(i,j) = u_opt(1)/sqrt(u_opt(1)^2 + u_opt(2)^2);
                        self.uy(i,j) = u_opt(2)/sqrt(u_opt(1)^2 + u_opt(2)^2);
                        %for not-normalized vector field
                        %plotux(i,j) = u_opt(1);
                        %plotuy(i,j) = u_opt(2);
                    end
                end
            end
        end
        
        function make_plot(self)
            % PLOT ELLIPTICAL OBSTACLE
            L = chol(self.Q);
            ellipse_pts = linspace(0,2*pi,300);
            z = [cos(ellipse_pts); sin(ellipse_pts)]*self.r; % creates circle of correct radius
            %R = rotation_matrix(pi/4);
            ellipse = L^-1 * z + self.q';    % stretch and translate circle appropriately

            % PLOT GENERATED TRAJECTORY and VECTOR FIELD
            axis(self.limits);
            hold on
            plot(ellipse(1,:),ellipse(2,:),'black','LineWidth',3)
            quiver(self.ax,self.ay,self.ux,self.uy,'red')
            for i = 1:self.numtrajectories
                scatter(self.x(1,1,i),self.x(1,2,i),50,'blue','filled')
                plot(self.x(:,1,i),self.x(:,2,i),'blue','LineWidth',2)
            end
            %legend('Obstacle Boundary','Input Vector','Initial State','Trajectory')
        end
    end
end
