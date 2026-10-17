% =============================================================================
% newmark_nonlinear.m — Nonlinear Newmark Time Integration Method
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Solves the nonlinear equation of motion using a predictor-corrector
%   scheme based on the Newmark method with Newton-Raphson iterations:
%
%       M·ü + C·u̇ + K·u + Knl·u³ = F(t)
%
%   At each time step, displacement is predicted, then corrected iteratively
%   using the tangent stiffness matrix until the residual converges.
%
% Inputs:
%   M     - Mass matrix                          [Ne x Ne]
%   K     - Linear stiffness matrix              [Ne x Ne]
%   Knl   - Cubic stiffness per DOF              [scalar or Ne x 1]
%   C     - Damping matrix                       [Ne x Ne]
%   F     - External force matrix                [Ne x Np]
%   u0    - Initial displacement vector          [Ne x 1]
%   u0p   - Initial velocity vector              [Ne x 1]
%   dt    - Time step                            [scalar, s]
%   Np    - Number of time steps                 [integer]
%   gamma - Newmark parameter (>= 0.5)           [scalar]
%   beta  - Newmark parameter (>= 0.25*(0.5+gamma)^2) [scalar]
%   tol   - Convergence tolerance for residual   [scalar]
%
% Outputs:
%   u     - Displacement time history            [Ne x Np]
%   up    - Velocity time history                [Ne x Np]
%   upp   - Acceleration time history            [Ne x Np]
%
% Recommended parameters:
%   Average acceleration (unconditionally stable): gamma = 0.5, beta = 0.25
%   Linear acceleration:                           gamma = 0.5, beta = 1/6
%
% Reference:
%   Geradin & Rixen (2014). Mechanical Vibrations. Wiley. Ch. 9.
%   Gass, N.L. (2025). Nonlinear modeling of clunk noise in vehicle
%   transmission systems. MSc Dissertation, UNICAMP.
% =============================================================================

function [u, up, upp] = newmark_nonlinear(M, K, Knl, C, F, u0, u0p, dt, Np, gamma, beta, tol)

%% --- Initial Setup ---

Ne = size(u0, 1);                           % Number of degrees of freedom

% Initial acceleration from equation of motion: M*u0pp = F0 - C*u0p - P(u0)
P0   = C*u0p(:,1) + K*u0(:,1) + Knl.*(u0(:,1)).^3;
u0pp = M \ (F(:,1) - P0);

% Newmark integration constants
a1 = (1 - gamma) * dt;                      % Predictor coeff. for velocity
a2 = (0.5 - beta) * dt^2;                   % Predictor coeff. for displacement
a3 = gamma / (beta * dt);                   % Corrector coeff. for velocity
a4 = 1 / (beta * dt^2);                     % Corrector coeff. for acceleration

%% --- Memory Allocation ---

u   = zeros(Ne, Np);
up  = zeros(Ne, Np);
upp = zeros(Ne, Np);

% Extended arrays used in the predictor-corrector loop
ue   = zeros(Ne, Np+1);   % Predicted displacement
upe  = zeros(Ne, Np+1);   % Predicted velocity
uppe = zeros(Ne, Np+1);   % Predicted acceleration
Pe   = zeros(Ne, Np+1);   % Internal force vector P(u, u̇)
Psi  = zeros(Ne, Np+1);   % Residual Ψ = M*ü + P - F
dtu  = zeros(Ne, Np+1);   % Newton-Raphson displacement correction δu

% Assign initial conditions
u(:,1)   = u0;
up(:,1)  = u0p;
upp(:,1) = u0pp;

% Initial residual (set > tol to enter while loop on first step)
Psi(:,1) = M*u0pp(:,1) + P0 - F(:,1);

%% --- Tangent Stiffness Function ---
% For cubic nonlinearity: Kt = K + diag(3*Knl.*u²)
% For linear system (Knl = 0): Kt = K
if any(Knl ~= 0)
    Kt = @(u_curr) K + diag(3 * Knl .* u_curr.^2);   % cubic springs act on each DOF
else
    Kt = @(u_curr) K;
end

%% --- Time Integration Loop ---

for i = 1:Np-1

    % ---- PREDICTOR STEP ----
    % Estimate state at t_{i+1} assuming zero acceleration increment
    uppe(:,i+1) = 0;
    upe(:,i+1)  = up(:,i)  + a1 * upp(:,i);
    ue(:,i+1)   = u(:,i)   + dt * up(:,i) + a2 * upp(:,i);

    % Initialize residual to enter corrector loop
    Psi(:,i+1) = tol + 1;

    % ---- CORRECTOR STEP (Newton-Raphson) ----
    iter = 0;
    while norm(Psi(:,i+1)) > tol
        iter = iter + 1;
        if iter > 100
            warning('newmark_nonlinear: Newton-Raphson did not converge at step %d (iter > 100).', i);
            break;
        end

        % Internal force at predicted state
        Pe(:,i+1) = C*upe(:,i+1) + K*ue(:,i+1) + Knl.*(ue(:,i+1)).^3;

        % Residual: Ψ = M*ü* + P(u*, u̇*) - F
        Psi(:,i+1) = M*uppe(:,i+1) + Pe(:,i+1) - F(:,i+1);

        % Effective stiffness: S = Kt + (γ/βΔt)·C + (1/βΔt²)·M
        S = Kt(ue(:,i+1)) + a3*C + a4*M;

        % Newton-Raphson correction
        dtu(:,i+1) = -S \ Psi(:,i+1);

        % Update state with correction
        u(:,i+1)   = ue(:,i+1)   + dtu(:,i+1);
        up(:,i+1)  = upe(:,i+1)  + a3 * dtu(:,i+1);
        upp(:,i+1) = uppe(:,i+1) + a4 * dtu(:,i+1);

        % Recompute internal force and residual with corrected state
        Pe(:,i+1)  = C*up(:,i+1) + K*u(:,i+1) + Knl.*(u(:,i+1)).^3;
        Psi(:,i+1) = M*upp(:,i+1) + Pe(:,i+1) - F(:,i+1);

        % Update predicted state for next iteration
        ue(:,i+1)   = u(:,i+1);
        upe(:,i+1)  = up(:,i+1);
        uppe(:,i+1) = upp(:,i+1);
    end

end

end
