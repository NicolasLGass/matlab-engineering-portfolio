% =============================================================================
% runge_kutta.m — 4th-Order Runge-Kutta Solver (RK4)
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Classic 4th-order Runge-Kutta (RK4) solver for first-order ODEs.
%   Used as an independent reference/validation for the Newmark methods.
%
%   Solves: dq/dt = func(t, q, F_input)
%
%   For a 2nd-order SDOF system M*x'' + C*x' + K*x = F(t),
%   define state vector q = [x; x'] and:
%     func = @(t,q,F) [q(2); (F - C*q(2) - K*q(1)) / M]
%
% Inputs:
%   func  - Function handle: func(t, q, u) -> dq/dt    [function handle]
%   t     - Time vector                                [1 x Np]
%   q0    - Initial state vector                       [Nq x 1]
%   F     - Scalar force input at each time step       [1 x Np]
%   h     - Time step size (must match t spacing)      [scalar, s]
%
% Outputs:
%   qf    - State time history, rows = time steps      [Np x Nq]
%
% Usage example (SDOF):
%   zdot = @(t,z,F) [z(2); (F - C*z(2) - K*z(1)) / M];
%   qf   = runge_kutta(zdot, t, [x0; v0], F_vec, dt);
%   x    = qf(:,1)';   % Displacement
%   v    = qf(:,2)';   % Velocity
%
% Reference:
%   Press et al. (2007). Numerical Recipes: The Art of Scientific Computing.
%   Cambridge University Press. Section 17.1.
% =============================================================================

function [qf] = runge_kutta(func, t, q0, F, h)

Np = length(t);
Nq = length(q0);

qf    = zeros(Np, Nq);
qf(1,:) = q0';
q     = q0;

for i = 1:Np-1

    % Classic RK4 coefficients
    k1 = func(t(i),        q,                 F(i));
    k2 = func(t(i) + 0.5*h, q + 0.5*k1*h,    F(i));
    k3 = func(t(i) + 0.5*h, q + 0.5*k2*h,    F(i));
    k4 = func(t(i) + h,     q + k3*h,         F(i));

    % Weighted average update
    q = q + ((k1 + 2*k2 + 2*k3 + k4) / 6) * h;
    qf(i+1,:) = q';

end

end
