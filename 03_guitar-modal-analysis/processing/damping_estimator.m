% =============================================================================
% damping_estimator.m — Logarithmic Decrement Damping Estimation
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          MEAV — UFSM (2019)
%
% Description:
%   Estimates the viscous damping ratio (zeta) from a free-decay vibration
%   signal using the logarithmic decrement method.
%
%   Method:
%     1. Extract decay envelope via find_decay()
%     2. Detect the free-decay peaks (one per damped period)
%     3. wd from peak spacing; zeta from a least-squares fit of ln(peaks) vs time
%
%   Logarithmic decrement:
%     delta = ln(x_n / x_{n+p}) / p
%     zeta  = delta / sqrt(4*pi^2 + delta^2)  ≈ delta / (2*pi)  for small zeta
%
% Inputs (workspace):
%   accz  - Acceleration signal (free decay) [V or m/s²]
%   Fs    - Sampling rate [Hz]
%
% Outputs (workspace):
%   zeta_est  - Estimated damping ratio [-]
%   wn_est    - Estimated natural frequency [rad/s]
%   fn_est    - Estimated natural frequency [Hz]
%
% =============================================================================

%% =========================================================================
%  SECTION 1 — LOAD DATA
%  Load your measured data — must contain accz and Fs
% ==========================================================================

% load(fullfile('data', 'free_decay.mat'));   % your file: accz (free-decay signal), Fs

% Verify
if ~exist('accz', 'var') || ~exist('Fs', 'var')
    error('Load your data file first. Required variables: accz, Fs');
end

%% =========================================================================
%  SECTION 2 — EXTRACT DECAY ENVELOPE
% ==========================================================================

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));   % find_decay, local_peaks

n_smooth = 200;     % Envelope smoothing: larger = smoother (tune per signal)

[env_upper, sig_clean] = find_decay(accz, n_smooth);

t_sig = (0 : length(sig_clean)-1) / Fs;
t_env = (0 : length(env_upper)-1) / Fs;

%% =========================================================================
%  SECTION 3 — ESTIMATE NATURAL FREQUENCY FROM PEAKS
% ==========================================================================

[pks, locs] = local_peaks(sig_clean, 1, [], max(sig_clean)*0.1);   % = findpeaks(..., 'MinPeakProminence', 10% of max)

if length(locs) >= 2
    T_periods = diff(locs) / Fs;        % Period between consecutive peaks [s]
    Td_avg    = mean(T_periods);        % Average damped period [s]
    fd_est    = 1 / Td_avg;            % Damped natural frequency [Hz]
    wd_est    = 2*pi * fd_est;         % [rad/s]
else
    error('Not enough peaks found. Check signal quality or reduce MinPeakProminence.');
end

fprintf('Damped natural frequency: fd = %.3f Hz  (wd = %.3f rad/s)\n', fd_est, wd_est);

%% =========================================================================
%  SECTION 4 — ESTIMATE DAMPING RATIO (Logarithmic Decrement)
% ==========================================================================

% Peaks are one damped period apart, so ln(x_k) = ln(x_0) - (zeta*wn)*t_k.
% A least-squares line through ln(pks) gives the decay rate; times Td it is
% the log decrement averaged over all detected cycles.
coef     = polyfit(locs(:)/Fs, log(pks(:)), 1);   % slope = -zeta*wn
delta    = -coef(1) * Td_avg;                     % Log decrement per cycle
zeta_est = delta / sqrt(4*pi^2 + delta^2);

wn_est = wd_est / sqrt(1 - zeta_est^2);
fn_est = wn_est / (2*pi);

fprintf('Damping ratio (log. decrement): zeta = %.5f\n', zeta_est);
fprintf('Natural frequency: fn = %.3f Hz  (wn = %.3f rad/s)\n', fn_est, wn_est);

%% =========================================================================
%  SECTION 5 — FIT VISUALIZATION
% ==========================================================================

A0  = env_upper(1);
env_fit = A0 * exp(-zeta_est * wn_est * t_env);

figure('Name', 'Damping Estimation — Logarithmic Decrement');
subplot(2,1,1);
plot(t_sig, sig_clean, 'b', 'LineWidth', 0.8, 'DisplayName', 'Signal'); hold on;
plot(t_env, env_upper, 'r', 'LineWidth', 1.5, 'DisplayName', 'Envelope (peak)');
plot(t_env, env_fit,  'k--', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Fit: A_0 e^{-\\zeta\\omega_n t},  \\zeta=%.4f', zeta_est));
grid on; xlabel('Time [s]'); ylabel('Amplitude');
title('Free-Decay Signal and Envelope Fit');
legend('Location', 'northeast');

subplot(2,1,2);
semilogy(t_env, env_upper, 'r', 'LineWidth', 1.5, 'DisplayName', 'Envelope');
hold on;
semilogy(t_env, env_fit, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Exponential fit');
grid on; xlabel('Time [s]'); ylabel('Amplitude (log scale)');
title(sprintf('Envelope Decay  |  f_n = %.2f Hz, \\zeta = %.4f', fn_est, zeta_est));
legend;
