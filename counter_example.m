method = 'standard';
%method = 'jankovic';

% SYSTEM AND INTEGRATION PARAMETERS
sysdim = 2; % system dimension
numtrajectories = 2; % number of trajectories to simulate
tsteps = 3000; % number of timesteps to execute
stepsize = 1E-3; % step size between iterations
u = zeros(tsteps,sysdim,numtrajectories); % create vector of inputs
x = zeros(tsteps,sysdim,numtrajectories); % create vector of states
t = linspace(0,stepsize*tsteps,tsteps); % generate vector of times

% INITIAL STATES AND INPUTS
initial_states = [0 3.75; -3 3.5]';
initial_inputs = zeros(numtrajectories,sysdim)';
x(1,:,:) = initial_states;
u(1,:,:) = initial_inputs;

if method == 'standard'
    % GENERATE TRAJECTORY
    for j = 1:numtrajectories
        for i = 2:tsteps
            % SOLVE STANDARD CLF - CBF OPTIMIZATION PROBLEM
            u(i,:,j) = CLF_CBF_QP(x(i-1,:,j));

            % PERFORM INTEGRATION STEP
            % Using ode45 for laziness
            tspan = [t(i-1),t(i)]; % time interval for current iteration
            [sol_t, sol_x] = ode45(@(t,y) xdot(t,y,u(i,:,j)), tspan, x(i-1,:,j)); % use ode45 to integrate
            x(i,:,j) = sol_x(end,:);
        end
    end
elseif method == 'jankovic'
    % GENERATE TRAJECTORY
    for j = 1:numtrajectories
        for i = 2:tsteps
            % JANKOVIC's MODIFIED CLF - CBF OPTIMIZATION PROBLEM    
            u(i,:,j) = Jankovic_CLF_CBF_QP(x(i-1,:,j));

            % PERFORM INTEGRATION STEP
            % Using ode45 for laziness
            tspan = [t(i-1),t(i)]; % time interval for current iteration
            [sol_t, sol_x] = ode45(@(t,y) xdot(t,y,u(i,:,j)), tspan, x(i-1,:,j)); % use ode45 to integrate
            x(i,:,j) = sol_x(end,:);
        end
    end
end

% GENERATE VECTOR FIELD
vector_field;

% PLOT RESULTS
plotting;

% EQUATIONS OF MOTION
% dx/dt = f(x) + g(x)u = u
function odefcn = xdot(t,x,u)
[A, B] = get_linear_dynamics();
odefcn = A*x + B*(u');
end