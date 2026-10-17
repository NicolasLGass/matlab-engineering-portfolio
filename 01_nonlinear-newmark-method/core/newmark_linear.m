% =============================================================================
% newmark_linear.m — Linear Newmark Time Integration Method
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Solves the linear equation of motion:
%
%       M·ü + C·u̇ + K·u = F(t)
%
%   using the standard Newmark method with constant effective stiffness
%   matrix (no iterations required — O(N) per step).
%
% Inputs:
%   M     - Mass matrix                          [Ne x Ne]
%   K     - Stiffness matrix                     [Ne x Ne]
%   C     - Damping matrix                       [Ne x Ne]
%   F     - External force matrix                [Ne x Np]
%   u0    - Initial displacement                 [Ne x 1]
%   u0p   - Initial velocity                     [Ne x 1]
%   u0pp  - Initial acceleration                 [Ne x 1]
%   dt    - Time step                            [scalar, s]
%   Np    - Number of time steps                 [integer]
%   gamma - Newmark gamma parameter (>= 0.5)
%   beta  - Newmark beta  parameter (>= 0.25*(0.5+gamma)^2)
%
% Outputs:
%   u     - Displacement time history            [Ne x Np]
%   up    - Velocity time history                [Ne x Np]
%   upp   - Acceleration time history            [Ne x Np]
%
% Note:
%   For gamma = 0.5, beta = 0.25 (average acceleration): unconditionally stable.
%   For gamma = 0.5, beta = 1/6  (linear acceleration) : conditionally stable.
%
% Reference:
%   Geradin & Rixen (2014). Mechanical Vibrations. Wiley.
%   Bathe & Wilson (1976). Numerical Methods in Finite Element Analysis.
% =============================================================================

function [u, up, upp] = newmark_linear(M, K, C, F, u0, u0p, u0pp, dt, Np, gamma, beta)

%% --- Newmark Constants ---

a0 = 1 / (beta * dt^2);
a1 = gamma / (beta * dt);
a2 = 1 / (beta * dt);
a3 = 1/(2*beta) - 1;
a4 = gamma/beta - 1;
a5 = (dt/2) * (gamma/beta - 2);
a6 = dt * (1 - gamma);
a7 = gamma * dt;

% Effective (constant) stiffness matrix
Keff = K + a0*M + a1*C;

%% --- Memory Allocation ---

Ne  = size(u0, 1);
u   = zeros(Ne, Np);
up  = zeros(Ne, Np);
upp = zeros(Ne, Np);
Feff = zeros(Ne, Np);

% Initial conditions
u(:,1)   = u0;
up(:,1)  = u0p;
upp(:,1) = u0pp;

%% --- Time Integration Loop ---

for i = 1:Np-1

    % Effective force vector
    Feff(:,i+1) = F(:,i+1) ...
        + M * (a0*u(:,i) + a2*up(:,i) + a3*upp(:,i)) ...
        + C * (a1*u(:,i) + a4*up(:,i) + a5*upp(:,i));

    % Solve for displacement
    u(:,i+1) = Keff \ Feff(:,i+1);

    % Update acceleration and velocity
    upp(:,i+1) = a0*(u(:,i+1) - u(:,i)) - a2*up(:,i) - a3*upp(:,i);
    up(:,i+1)  = up(:,i) + a6*upp(:,i) + a7*upp(:,i+1);

end

end
