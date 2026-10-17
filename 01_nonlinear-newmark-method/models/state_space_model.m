% =============================================================================
% state_space_model.m — State-Space Formulation and lsim/ode Reference
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Builds the state-space representation of the linear SDOF oscillator
%   and computes the response using MATLAB's lsim() for comparison with
%   the Newmark and Runge-Kutta solvers.
%
%   State-space form: ẋ = A*x + B*u,   y = C*x + D*u
%   State vector: x = [displacement; velocity]
%
% Requires workspace variables: M, K, C, F, t, x0, v0, ccase
%   (defined in system_parameters.m and force_vector.m)
%
% Outputs (added to workspace):
%   r_ss  - Response from lsim [1 x Np], displacement [m]
% =============================================================================

%% --- State-Space Matrices ---

A_ss = [0,    1;
       -K/M, -C/M];

B_ss = [0;
        1/M];

C_ss = [1, 0];   % Output: displacement only

D_ss = 0;

if exist('OCTAVE_VERSION', 'builtin'), pkg load control; end   % Octave: control package
sys  = ss(A_ss, B_ss, C_ss, D_ss);
x_ini = [x0; v0];

%% --- Compute Response ---

if ccase == 1
    % Free vibration — initial conditions only, no input
    r_ss = initial(sys, x_ini, t)';
else
    % Forced response — uses force matrix F (row vector required by lsim)
    r_ss = lsim(sys, F(1,:)', t, x_ini)';
end

fprintf('State-space model computed (%d time steps).\n', length(t));
