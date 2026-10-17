% =============================================================================
% engine_mount_analysis.m — Rubber Engine Mount: Parameter Estimation
% =============================================================================
% Author : Nicolas Londero Gass — 16/08/2018
%          Mechanical Engineering — UFSM
%
% Description:
%   Estimates stiffness and transmissibility of cylindrical rubber engine
%   mounts (coxim) for vibration isolation of a 40 Hz excitation source.
%
%   Design criterion:
%     The mount resonance frequency must satisfy:
%       f_r = (1/2pi) * sqrt(K_total / M) < 40/sqrt(2) ≈ 28.28 Hz
%     This ensures at least -3 dB attenuation at the target frequency.
%
%   Stiffness model for a cylindrical rubber isolator under axial compression:
%     K_i = E*(1+j*eta)*A / L
%     K_total = 4 * K_i    (4 mounts in parallel)
%
%   where:
%     E    = Dynamic Young's modulus of rubber   [Pa]
%     eta  = Loss factor                          [-]
%     A    = Contact area = pi*D^2/4             [m²]
%     L    = Isolator height                      [m]
%
%   Force transmissibility:
%     T(omega) = |K_t / (K_t - M*omega^2)|
%
% Note on data:
%   Geometric and material parameters are from an academic estimation study.
%   Supplier, product and material trade names are intentionally omitted.
%
% Shore A to Young's modulus conversion (approximate):
%   Shore 40: E ≈ 1.58 MPa  |  Shore 50: E ≈ 2.24 MPa
%   Shore 60: E ≈ 3.55 MPa  |  Shore 70: E ≈ 5.62 MPa
% =============================================================================

clear; close all; clc; tic;

%% =========================================================================
%  SECTION 1 — SYSTEM PARAMETERS
% ==========================================================================

M   = 26.4;             % Supported mass              [kg]
eta = 0.1;              % Rubber loss factor          [-]

f_target = 40;          % Excitation frequency        [Hz]
f_r_max  = f_target / sqrt(2);   % Max resonance freq [Hz] (for -3 dB)

f   = 1:0.1:100;        % Frequency vector            [Hz]
w   = 2 * pi * f;       % Angular frequency           [rad/s]

%% =========================================================================
%  SECTION 2A — GEOMETRIC CONFIGURATIONS (fixed Shore hardness)
%   Exploring D and L for E = 10 MPa rubber
% ==========================================================================

D_geo = [0.024, 0.023, 0.0195, 0.0155, 0.015, 0.010, 0.0335];  % Diameter [m]
L_geo = [0.022, 0.020, 0.015,  0.010,  0.020, 0.020, 0.032 ];  % Length   [m]
E_geo = [0.01,  0.01,  0.01,   0.01,   0.01,  0.05,  0.0038] * 1e9; % Young's [Pa]

labels_geo = {
    'D=24mm, L=22mm, E=10MPa',
    'D=23mm, L=20mm, E=10MPa',
    'D=19.5mm, L=15mm, E=10MPa',
    'D=15.5mm, L=10mm, E=10MPa',
    'D=15mm, L=20mm, E=10MPa',
    'D=10mm, L=20mm, E=50MPa',
    'Reference mount: D=33.5mm, L=32mm, E=3.8MPa'
};

[Tf_geo, fr_geo] = compute_transmissibility(D_geo, L_geo, E_geo, eta, M, w);

% Print resonance frequencies
fprintf('--- Geometric configurations (resonance frequencies) ---\n');
for i = 1:length(D_geo)
    fprintf('  Config %d: f_r = %.2f Hz  (target < %.2f Hz)\n', i, fr_geo(i), f_r_max);
end

%% =========================================================================
%  SECTION 2B — SHORE HARDNESS SWEEP (fixed geometry)
%   D=25mm, L=21mm, varying Shore A 40/50/60/70
% ==========================================================================

D_sh  = repmat(0.025, 1, 4);   % Diameter [m]
L_sh  = repmat(0.021, 1, 4);   % Length   [m]
E_sh  = [1.5848, 2.2387, 3.5481, 5.6234] * 1e6;  % Young's moduli [Pa]

labels_sh = {'Shore 40 (E=1.58 MPa)', 'Shore 50 (E=2.24 MPa)', ...
             'Shore 60 (E=3.55 MPa)', 'Shore 70 (E=5.62 MPa)'};

[Tf_sh, fr_sh] = compute_transmissibility(D_sh, L_sh, E_sh, eta, M, w);

fprintf('\n--- Shore hardness sweep ---\n');
for i = 1:4
    fprintf('  %s: f_r = %.2f Hz\n', labels_sh{i}, fr_sh(i));
end

%% =========================================================================
%  SECTION 3 — PLOTS
% ==========================================================================

figure('Name', 'Engine Mount — Geometric Configurations', 'Position', [50 100 900 500]);
semilogy(f, Tf_geo, 'LineWidth', 1.6);
hold on; grid on;
% Reference lines (portable replacement for xline/yline)
plot([f_target f_target], [0.01 100], 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot([f_r_max f_r_max],   [0.01 100], 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot([1 100], [1 1], 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
text(f_target, 30, ' 40 Hz target'); text(f_r_max, 30, ' 28.3 Hz limit', 'Color', 'r');
xlabel('Frequency [Hz]'); ylabel('Transmissibility [-]');
title('Force Transmissibility — Geometric Configurations');
legend(labels_geo, 'Location', 'southeast', 'FontSize', 8);
xlim([1 100]); ylim([0.01 100]);

figure('Name', 'Engine Mount — Shore Hardness', 'Position', [50 150 900 450]);
semilogy(f, Tf_sh, 'LineWidth', 1.8);
hold on; grid on;
% Reference lines (portable replacement for xline/yline)
plot([f_target f_target], [0.01 100], 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot([f_r_max f_r_max],   [0.01 100], 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
plot([1 100], [1 1], 'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
text(f_target, 30, ' 40 Hz target'); text(f_r_max, 30, ' 28.3 Hz limit', 'Color', 'r');
xlabel('Frequency [Hz]'); ylabel('Transmissibility [-]');
title('Force Transmissibility — Shore Hardness Sweep (D=25mm, L=21mm)');
legend(labels_sh, 'Location', 'southeast');
xlim([1 100]); ylim([0.01 100]);

toc;
