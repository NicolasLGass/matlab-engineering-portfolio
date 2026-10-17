% =============================================================================
% validate_amesim.m — Nonlinear Newmark vs. Siemens Simcenter AMESim
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Reproduces the comparison between newmark_nonlinear.m and the same
%   SDOF system modelled in Simcenter AMESim (mass + damper + cubic spring):
%
%       m·ü + c·u̇ + k·u + Knl·u³ = F0·sin(2π·f·t)
%
%   m = 20 kg, k = 20 N/m, c = 4 N·s/m, f = 0.2 Hz, u(0) = u̇(0) = 0
%
%   Four cases, from linear to strongly nonlinear (hardening spring):
%       Knl [N/m³] :   0     10    100    500
%       F0  [N]    :  10     10    100   1000
%
%   The AMESim results (data/amesim_harmonic_0p2Hz.csv) were exported at
%   Δt = 1 ms and are stored decimated to 10 ms to keep the file small.
%   Newmark runs at Δt = 1 ms and is compared at the AMESim sample times.
%
% Output:
%   Table with peak and RMS errors (normalized by the AMESim response) and
%   a figure with the time histories and error for each case.
% =============================================================================

clear; close all; clc;

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'core'));

%% --- System and integration parameters ---
m = 20;  k = 20;  c = 4;            % [kg], [N/m], [N·s/m]
f_exc = 0.2;                        % Excitation frequency [Hz]

dt    = 1e-3;                       % Newmark time step [s]
Tf    = 20;                         % Duration [s]
t     = 0:dt:Tf;
Np    = length(t);
gamma = 0.5;  beta = 0.25;          % Average acceleration
tol   = 1e-10;                      % Newton-Raphson tolerance

cases = struct('Knl', {0, 10, 100, 500}, 'F0', {10, 10, 100, 1000});

%% --- AMESim reference ---
ref   = dlmread(fullfile(here, 'data', 'amesim_harmonic_0p2Hz.csv'), ',', 1, 0);  % skip header (MATLAB/Octave)
t_as  = ref(:,1);
idx   = round(t_as/dt) + 1;         % Newmark samples matching AMESim times

%% --- Run and compare ---
fprintf('  Knl [N/m³]   F0 [N]   max|x| [m]   peak error   RMS error\n');
fprintf('  ----------------------------------------------------------\n');

figure('Name', 'Newmark NL vs AMESim', 'Position', [80 80 1000 700]);
for n = 1:numel(cases)
    Knl = cases(n).Knl;
    F   = cases(n).F0 * sin(2*pi*f_exc*t);

    u   = newmark_nonlinear(m, k, Knl, c, F, 0, 0, dt, Np, gamma, beta, tol);
    x_nm = u(idx).';
    x_as = ref(:, n+1);

    err      = x_nm - x_as;
    peak_err = max(abs(err)) / max(abs(x_as)) * 100;
    rms_err  = sqrt(mean(err.^2) / mean(x_as.^2)) * 100;
    fprintf('  %8g   %7g   %10.4f   %9.4f %%  %8.4f %%\n', ...
            Knl, cases(n).F0, max(abs(x_as)), peak_err, rms_err);

    subplot(numel(cases), 2, 2*n-1);
    plot(t_as, x_as, 'k', 'LineWidth', 2.2); hold on;
    plot(t_as, x_nm, '--', 'Color', [0 0.447 0.741], 'LineWidth', 1.4);
    grid on; ylabel('x [m]');
    title(sprintf('K_{nl} = %g N/m^3, F_0 = %g N', Knl, cases(n).F0));
    if n == 1, legend('AMESim', 'Newmark NL', 'Location', 'northeast'); end
    if n == numel(cases), xlabel('Time [s]'); end

    subplot(numel(cases), 2, 2*n);
    plot(t_as, err, 'Color', [0.85 0.325 0.098]); grid on;
    ylabel('x_{NM} - x_{AS} [m]');
    title(sprintf('Peak error %.4f %%', peak_err));
    if n == numel(cases), xlabel('Time [s]'); end
end
if exist('sgtitle') ~= 0                    % not available in older Octave versions
    sgtitle('Nonlinear Newmark (\Deltat = 1 ms) vs Simcenter AMESim');
end
