% =============================================================================
% main_impact_pair.m — Example: Nonlinear Impact Pair (Gear Backlash)
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Self-contained example demonstrating the nonlinear Newmark solver for
%   a single-DOF impact pair model (gear rattle / clunk noise).
%
%   The impact pair adds a cubic nonlinear stiffness term:
%       m*x'' + c*x' + k*x + Knl*x^3 = F(t)
%
%   Three solvers are compared:
%     1. Nonlinear Newmark (newmark_nonlinear.m)
%     2. Linear Newmark    (newmark_linear.m) — baseline
%     3. Runge-Kutta 4     (runge_kutta.m)   — reference
%
%   The three Knl values (0, 1, 10) reproduce the parametric study from
%   Chapter 4 of the MSc dissertation.
%
% Reference:
%   Gass, N.L. (2025). Nonlinear modeling of clunk noise in vehicle
%   transmission systems. MSc Dissertation, UNICAMP.
%   Comparin & Singh (1988). J. Sound and Vibration, 127(2), 337–358.
% =============================================================================

clear; close all; clc;

%% --- Add solver paths ---
here = fileparts(mfilename('fullpath'));            % works from any current folder
addpath(fullfile(here, '..', 'core'), fullfile(here, '..', 'models'));

%% =========================================================================
%  SECTION 1 — SYSTEM PARAMETERS
% ==========================================================================

% Physical parameters
m1 = 20;        % Mass          [kg]
k1 = 20;        % Stiffness     [N/m]
c1 = 4;         % Damping       [N·s/m]

M = m1;
K = k1;
C = c1;
Ne = 1;

% Natural frequency
wn   = sqrt(K/M);
fn   = wn / (2*pi);
fprintf('Natural frequency: %.3f Hz\n', fn);

%% =========================================================================
%  SECTION 2 — SIMULATION PARAMETERS
% ==========================================================================

dt    = 0.001;              % Time step              [s]
Tf    = 40;                 % Simulation duration    [s]
t     = 0:dt:Tf;
Np    = length(t);

gamma = 0.5;                % Newmark gamma (average acceleration)
beta  = 0.25;               % Newmark beta
tol   = 1e-6;               % Newton-Raphson tolerance

%% =========================================================================
%  SECTION 3 — EXCITATION FORCE
% ==========================================================================
% Choose force type: ccase = 1 - Free vibration (no external force)
%                    ccase = 2 - Unit impulse approximation
%                    ccase = 3 - Constant step force
%                    ccase = 4 - Harmonic force: F = Fmag*sin(2*pi*exfreq*t)
%                    ccase = 5 - Step force ending at half-time

ccase  = 3;         % Constant step force
Fmag   = 10;        % Force magnitude [N]
timp   = 0;
timpf  = 0;
exfreq = fn;

[F, ~] = force_vector(t, ccase, Np, exfreq, Fmag, timp, timpf, dt, Ne);

%% =========================================================================
%  SECTION 4 — SOLVE FOR MULTIPLE Knl VALUES
% ==========================================================================

Knl_cases = [0, 1, 10];     % Nonlinear stiffness values to compare
x0 = 0;                     % Initial displacement [m]
v0 = 0;                     % Initial velocity     [m/s]

% Initial acceleration
a0 = (F(1,1) - C*v0 - K*x0) / M;

results = struct();

for idx = 1:length(Knl_cases)
    Knl = Knl_cases(idx);

    % Nonlinear Newmark
    [u_nl, ~, a_nl] = newmark_nonlinear(M, K, Knl, C, F, x0, v0, dt, Np, gamma, beta, tol);

    % Linear Newmark (Knl = 0 baseline)
    [u_lin, ~, a_lin] = newmark_linear(M, K, C, F, x0, v0, a0, dt, Np, gamma, beta);

    % Runge-Kutta (linear only — for reference)
    zdot = @(tt, z, Ft) [z(2); (Ft - C*z(2) - K*z(1)) / M];
    qf   = runge_kutta(zdot, t, [x0; v0], F(1,:), dt);

    results(idx).Knl   = Knl;
    results(idx).u_nl  = u_nl;
    results(idx).u_lin = u_lin;
    results(idx).u_rk  = qf(:,1)';
    results(idx).a_nl  = a_nl;
    results(idx).a_lin = a_lin;
end

%% =========================================================================
%  SECTION 5 — PLOTS
% ==========================================================================

colors = {'#0072BD', '#D95319', '#77AC30'};
t_plot = t;

figure('Name', 'Displacement — Newmark NL vs Linear vs RK4', ...
       'Position', [100 100 900 600]);

for idx = 1:length(Knl_cases)
    subplot(length(Knl_cases), 1, idx);
    plot(t_plot, results(idx).u_lin, 'k--', 'LineWidth', 1, 'DisplayName', 'Linear Newmark');
    hold on; grid on;
    plot(t_plot, results(idx).u_rk,  '-.', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'RK4 (linear)');
    plot(t_plot, results(idx).u_nl,  '-',  'Color', colors{1}, 'LineWidth', 1.8, 'DisplayName', sprintf('Nonlinear Newmark (Knl=%.0f)', results(idx).Knl));
    xlabel('Time [s]');
    ylabel('Displacement [m]');
    title(sprintf('Response — K_{nl} = %.0f N/m³', results(idx).Knl));
    legend('Location', 'northeast');
    xlim([0, Tf]);
end
if exist('sgtitle') ~= 0, sgtitle('Impact Pair: Nonlinear vs Linear Response', 'FontSize', 13, 'FontWeight', 'bold'); end   % not in older Octave

figure('Name', 'Acceleration comparison', 'Position', [100 150 900 300]);
plot(t_plot, results(2).a_lin, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Linear Newmark');
hold on; grid on;
plot(t_plot, results(2).a_nl,  '-.', 'Color', colors{1}, 'LineWidth', 1.8, 'DisplayName', 'Nonlinear Newmark (Knl=1)');
xlabel('Time [s]'); ylabel('Acceleration [m/s²]');
title('Acceleration Response (Knl = 1 N/m³)');
legend; xlim([0, Tf]);
