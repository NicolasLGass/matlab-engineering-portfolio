% =============================================================================
% data_acquisition.m — Guitar Impact Measurement (NI cDAQ)
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          MEAV — UFSM (2019)
%
% Description:
%   Single-shot acquisition of a tri-axial accelerometer and an impact
%   force transducer on an acoustic guitar body. Computes the H1 FRF and
%   coherence for each axis and saves the time signals for the free-decay
%   damping analysis (processing/damping_estimator.m).
%
% Hardware:
%   - NI cDAQ-9234, 4 channels, 24-bit
%   - 3 × IEPE accelerometers (X, Y, Z)
%   - 1 × IEPE force transducer (impact hammer)
%
% Requirements:
%   - MATLAB Data Acquisition Toolbox, Signal Processing Toolbox
%   - NI-DAQmx driver
%
% Output:
%   <SAVE_FILE>.mat with accx, accy, accz, force, t, Fs, H, Co, freqh
% =============================================================================

clear; format long; clc; close all;
if exist('OCTAVE_VERSION', 'builtin')
    error('This acquisition script needs MATLAB + Data Acquisition Toolbox (NI-DAQmx); it does not run in Octave.');
end

%% =========================================================================
%  SECTION 1 — USER SETTINGS (edit here)
% ==========================================================================

SAVE_FILE = fullfile(pwd, 'data', 'guitar_impact_01.mat');   % Output file
Fs    = 51200;                  % Sampling rate [Hz]
T_rec = 5;                      % Record duration [s]
NFFT  = 1 * Fs;                 % Welch segment length (1 s)

%% =========================================================================
%  SECTION 2 — TRANSDUCER SENSITIVITIES
% ==========================================================================

sens_accx  = 1.029e-3;          % [V/(m/s²)] — channel 0
sens_accy  = 1.037e-3;          % [V/(m/s²)] — channel 1
sens_accz  = 1.050e-3;          % [V/(m/s²)] — channel 2
sens_force = 20.35e-3;          % [V/N]      — channel 3

cal = 1 ./ [sens_accx, sens_accy, sens_accz, sens_force];

%% =========================================================================
%  SECTION 3 — ACQUISITION
% ==========================================================================

s = daq.createSession('ni');
s.Rate              = Fs;
s.DurationInSeconds = T_rec;
addAnalogInputChannel(s, 'cDAQ1Mod1', [0, 1, 2, 3], 'IEPE');

disp('Start of recording — strike the guitar now.');
[raw, t] = s.startForeground();
disp('End of recording.');

accx  = raw(:,1) * cal(1);
accy  = raw(:,2) * cal(2);
accz  = raw(:,3) * cal(3);
force = raw(:,4) * cal(4);

%% =========================================================================
%  SECTION 4 — FRF (H1) AND COHERENCE
% ==========================================================================

win = hann(NFFT);
[H,  freqh] = tfestimate(force, [accx accy accz], win, NFFT/2, NFFT, Fs);
[Co, ~]     = mscohere(force, [accx accy accz], win, NFFT/2, NFFT, Fs);

%% =========================================================================
%  SECTION 5 — PLOTS
% ==========================================================================

figure('Name', 'Time signals');
subplot(2,1,1); plot(t, [accx accy accz]); grid on;
xlabel('Time [s]'); ylabel('a [m/s²]'); legend('acc x', 'acc y', 'acc z');
subplot(2,1,2); plot(t, force, 'k'); grid on;
xlabel('Time [s]'); ylabel('Force [N]');

axes_lbl = {'x', 'y', 'z'};
figure('Name', 'FRF and coherence');
for k = 1:3
    subplot(2,3,k);   loglog(freqh, abs(H(:,k))); grid on; xlim([20 2000]);
    xlabel('f [Hz]'); ylabel('|H| [m/s²/N]'); title(['H_' axes_lbl{k} '(f)']);
    subplot(2,3,k+3); semilogx(freqh, Co(:,k));   grid on; xlim([20 2000]);
    xlabel('f [Hz]'); ylabel('Coherence [-]');
end

%% =========================================================================
%  SECTION 6 — SAVE
% ==========================================================================

if ~exist(fileparts(SAVE_FILE), 'dir'), mkdir(fileparts(SAVE_FILE)); end
save(SAVE_FILE, 'accx', 'accy', 'accz', 'force', 't', 'Fs', 'H', 'Co', 'freqh');
fprintf('Data saved to: %s\n', SAVE_FILE);
