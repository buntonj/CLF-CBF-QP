README FILE

FILES CONTAINED IN THIS DIRECTORY:

- counter_example.m
- CLF_CBF_QP.m
- Jankovic_CLF_CBF_QP.m
- vector_field.m
- plotting.m

Details for each:


-counter_example.m

This matlab script will create, with one run, the trajectories for any number of (pre-specified) using either the standard CLF-CBF QP approach, or Jankovic's modified version depending on the value of the string 'method'. (either 'standard' or 'jankovic')

It then calls the function for generating the vector field. (vector_field.m)

To do so, it generates the obstacle parameters (Q, q, and r, which can be modified at the top of the script) and generates matrices to hold trajectories and inputs at each time step, based on a provided number of timesteps requested (parameter 'tsteps').  The discretization steps are also left as a parameter ('stepsize').

The set of initial states for each trajectory desired to be run are left to be input ('initial_states').  The number of trajectories simulated can be altered with this parameter, combined with the parameter 'numtrajectories'.  Initial inputs are assumed to be zero.

The script then iteratively solves the differential equation (whose equations of motion are described in the function 'xdot', stored at the bottom of the script and able to be altered if desired) using ode45 to numerically integrate over each interval.

After the trajectories are simulated, the vector fields are calculated with a call to 'vector_field.m', and all results are plotted with the script 'plotting.m'.


-CLF_CBF_QP.m

This function, given the current state and the parameters Q,q,r (defined in counter_example.m) for an elliptical obstacle, will compute the optimal input according to the standard CLF-CBF quadratic program approach.  It solves this by expressing the problem as a linearly constrained quadratic program and solving it with CPLEX, which is available with an academic license.  If a different solver is desired (say CVX, for example), the problem in question is constructed as:

 minimize 0.5*x'*H*x
 s. t.    Aineq*x <= bineq

and can be input to any QP solver of choice using these variables.
The function returns the optimal input if only one return variable is given (e.g. x = CLF_CBF_QP(x,Q,q,r)), but has the option to also return corresonding dual variable values and the value of the relaxation variable delta by instead writing:
[u, lambda, delta] = CLF_CBF_QP(x,Q,q,r).

In this function file, we find the functions for the safety and control lyapunov functions along with their Lie derivatives, which may be altered if desired.
In particular:

-- get_lyapunov_P()
	This helper function defines the positive definite matrix P for all other lie derivatives of the lyapunov function.  If a modification of V(x) is desired in only this manner, express it here.

-- V(x)
	This computes the value of the lyapunov function at a given point.  By default it computes it as x*P*x', with x a row vector and P being pulled from the helper function defined above.

-- LfV(x)
	Computes the lie derivative of V at state x along f.  Will need modified if V(x) or f is changed.

-- LgV(x)
	computes the lie derivative of V at state x along g.  Will need modified if V(x) or g is changed.

-- h(x,Q,q,r)
	Computes the safety function h at state x given elliptical safety border defined with ellipse parameters Q,q, and r.

-- Lfh(x,Q,q,r)
	Computes lie derivative of safety function h with elliptical parameters Q,q, and r along f.  Will need to be modified if h or f change.

-- Lgh(X,q,q,r)
	Computes lie derivaitve of safety function h with elliptical parameters Q,q, and r along g.  Will need to be modified if h or g change.

--gain_alpha and gain_gamma(x)
	Scalar K-infinity functions defining the gain on the constraint coefficients in the quadratic program.


-Jankovic_CLF_CBF_QP.m

Performs the same steps outlined in CLF_CBF_QP but with Jankovic's modifications, and does not contain gain_gamma and gain_alpha, instead holds jankovic_gamma(x), which computes jankovic's sign-dependent gain in the constraints.  All lie derivative functions here are independent of those defined in CLF_CBF_QP and will need to be modified independently. (I'll work on integrating these.)


- vector_field.m

Uses the functions CLF_CBF_QP.m or Jankovic_CLF_CBF_QP.m (depending on the set value of the parameter 'method') to compute the vector fields resulting from an elliptical obstacle with parameters defined by Q,q, and r as created in the file counter_example.m.

the x and y ranges and their intervals may be modified with the corresponding variables at the beginning of the file.


- plotting.m

Plots the results of vector fields and trajectories on a single image.  By default, the plot limits are set to the x and y ranges defined in vector_field.m.


This function performs the same task as CLF_CBF_QP.m but 
