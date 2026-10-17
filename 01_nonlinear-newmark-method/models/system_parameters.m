% =============================================================================
% system_parameters.m — System and Simulation Parameters
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%
% Description:
%   Defines all physical and numerical parameters for the 1-DOF nonlinear
%   oscillator example (Impact Pair / gear backlash model).
%   Consolidates the former Param_dt*.m and Properties_1dof*.m scripts.
%
%   Physical model:
%       m * x'' + c * x' + k * x + Knl * x^3 = F(t)
%
%   Run this script before calling the solvers. Adjust dt and Tf as needed
%   to explore accuracy vs. computation time trade-offs.
% =============================================================================

%% --- Physical Parameters ---

m1 = 20;        % Mass          [kg]
k1 = 20;        % Stiffness     [N/m]
c1 = 4;         % Damping       [N·s/m]

% System matrices (1-DOF)
M = m1;         % Mass matrix   [1 x 1]
K = k1;         % Stiffness     [1 x 1]
C = c1;         % Damping       [1 x 1]

Ne = size(M, 1);    % Number of DOFs

%% --- Derived Modal Parameters ---

wn   = sqrt(K/M);           % Natural frequency     [rad/s]
fn   = wn / (2*pi);         % Natural frequency     [Hz]
cc   = 2*sqrt(K*M);         % Critical damping      [N·s/m]
zeta = C / cc;              % Damping ratio         [-]
wd   = wn * sqrt(1-zeta^2); % Damped nat. freq.     [rad/s]

fprintf('--- System properties ---\n');
fprintf('  Natural frequency : %.4f rad/s  (%.4f Hz)\n', wn, fn);
fprintf('  Damping ratio     : %.4f\n', zeta);
fprintf('  Damped freq.      : %.4f rad/s\n', wd);

%% --- Simulation Parameters ---
% Adjust dt and Tf to balance accuracy and speed.
% Validation table (Impact Pair, Knl = 1, Fmag = 10):
%   dt = 0.01   -> disp error ~10.6%,  accel error ~282%
%   dt = 0.001  -> disp error ~0.107%, accel error ~1.51%
%   dt = 0.0001 -> disp error ~0.001%, accel error ~0.022%

dt  = 0.001;            % Time step          [s]
Tf  = 40;               % Final time         [s]
t   = 0:dt:Tf;          % Time vector        [s]
Np  = length(t);        % Number of steps

%% --- Newmark Integration Parameters ---
% Average acceleration: gamma=0.5, beta=0.25 (unconditionally stable)

gamma = 0.5;
beta  = 0.25;
tol   = 1e-6;           % Newton-Raphson convergence tolerance

%% --- Excitation Frequency ---

exfreq = fn;            % Excitation at natural frequency [Hz] (resonance)

%% --- Initial Conditions ---

x0 = 0;     % Initial displacement [m]
v0 = 0;     % Initial velocity     [m/s]
