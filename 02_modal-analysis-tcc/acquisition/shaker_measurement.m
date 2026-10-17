% =============================================================================
% shaker_measurement.m — Shaker (Random Excitation) Modal Measurement (NI cDAQ)
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          Undergraduate Thesis (TCC) — UFSM
%
% Description:
%   Data acquisition routine for experimental modal analysis using an
%   electrodynamic shaker driven by band-limited white noise. Three IEPE
%   accelerometers and one force transducer are read by an NI cDAQ-9234.
%
%   Acquisition strategy:
%     - White-noise excitation (same realization for every measurement),
%       played through the sound card into the shaker amplifier
%     - 6 excitation points × 2 accelerometer groups (3 channels each)
%     - H1 FRF and coherence by Welch averaging (Hann window, 50% overlap)
%     - Each record is accepted/rejected interactively
%
%   Signal chain:
%     oscillator() → sound card → amplifier → shaker → structure
%     NI cDAQ-9234 (IEPE) → MATLAB DAQ Toolbox → FRF (H1) → Saved .mat
%
% Hardware:
%   - NI cDAQ-9234, 4 channels, 24-bit
%   - 3 × IEPE accelerometers
%   - 1 × IEPE force transducer (shaker stinger)
%   - Electrodynamic shaker + power amplifier
%
% Requirements:
%   - MATLAB Data Acquisition Toolbox, Signal Processing Toolbox
%   - NI-DAQmx driver
%   - ../utils/oscillator.m (third-party signal generator)
%
% Usage:
%   1. Connect hardware and verify channel numbering.
%   2. Set SAVE_PATH and the settings below.
%   3. Run the script — follow on-screen prompts per measurement point.
% =============================================================================

clear; format long; clc; close all;
if exist('OCTAVE_VERSION', 'builtin')
    error('This acquisition script needs MATLAB + Data Acquisition Toolbox (NI-DAQmx); it does not run in Octave.');
end
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));

%% =========================================================================
%  SECTION 1 — USER SETTINGS (edit here)
% ==========================================================================

SAVE_PATH = fullfile(pwd, 'data', 'shaker');    % Output directory
Fs        = 5120;                               % Sampling rate [Hz]
%   Available rates for NI-9234: (131072/256) / n, n = 1..31

T_rec   = 40;                                   % Record duration [s]
N_EXC   = 6;                                    % Number of excitation points
ACC_GROUPS = 1:3:6;                             % Accelerometer groups (1-3, 4-6)
AMP     = 1;                                    % Excitation amplitude [0-1]

NFFT    = 1 * Fs;                               % Welch segment length (1 s)
N_OVLP  = NFFT/2;                               % 50% overlap
N_DFT   = 20 * NFFT;                            % Zero-padded DFT length (0.05 Hz)

%% =========================================================================
%  SECTION 2 — TRANSDUCER SENSITIVITIES
% ==========================================================================

sens_acc1  = 0.9964e-3;     % [V/(m/s²)] — Accelerometer channel 0
sens_acc2  = 1.036e-3;      % [V/(m/s²)] — Accelerometer channel 1
sens_acc3  = 1.025e-3;      % [V/(m/s²)] — Accelerometer channel 2
sens_force = 20.35e-3;      % [V/N]      — Force transducer channel 3

cal = 1 ./ [sens_acc1, sens_acc2, sens_acc3, sens_force];

%% =========================================================================
%  SECTION 3 — DAQ SESSION AND EXCITATION SIGNAL
% ==========================================================================

if ~exist(SAVE_PATH, 'dir'), mkdir(SAVE_PATH); end

s = daq.createSession('ni');
s.Rate              = Fs;
s.DurationInSeconds = T_rec;
addAnalogInputChannel(s, 'cDAQ1Mod1', [0, 1, 2, 3], 'IEPE');

% White noise with 0.1 s raised-cosine ramps, 1 s longer than the record
xs = oscillator('White Noise', T_rec + 1, [], 0.1, 0, Fs);
xs = AMP .* xs ./ max(abs(xs));

fprintf('DAQ initialized: Fs = %d Hz, T_rec = %d s, NFFT = %d\n', Fs, T_rec, NFFT);

Data = struct();
Data.SampRate      = Fs;
Data.MeasuringTime = T_rec;
Data.Sensibility   = struct('Acc1', sens_acc1, 'Acc2', sens_acc2, ...
                            'Acc3', sens_acc3, 'Force', sens_force);

%% =========================================================================
%  SECTION 4 — MEASUREMENT LOOP
% ==========================================================================

for IP = 1:N_EXC                    % Excitation point
    for I = ACC_GROUPS              % Accelerometer group (channels I..I+2)

        accepted = false;
        while ~accepted

            fprintf('\nAccelerometers at points %d, %d, %d | shaker at point %d\n', ...
                    I, I+1, I+2, IP);
            disp('Position the accelerometers and press any key to start...'); pause;

            % Play excitation and record simultaneously
            sound(xs, Fs);
            [raw, t] = s.startForeground();

            acc   = raw(:,1:3) .* cal(1:3);
            force = raw(:,4)   .* cal(4);

            % FRF estimation (H1) and coherence
            win = hann(NFFT);
            [H,  freqh] = tfestimate(force, acc, win, N_OVLP, N_DFT, Fs);
            [Co, ~]     = mscohere(force, acc, win, N_OVLP, N_DFT, Fs);

            % Quick-check plots
            figure(1); clf;
            subplot(3,1,1); plot(t, acc); grid on;
            xlabel('t [s]'); ylabel('a [m/s²]'); legend('Acc1','Acc2','Acc3');
            subplot(3,1,2); loglog(freqh, abs(H)); grid on; xlim([20 2000]);
            xlabel('f [Hz]'); ylabel('|H_1| [m/s²/N]');
            subplot(3,1,3); semilogx(freqh, Co); grid on; xlim([20 2000]);
            xlabel('f [Hz]'); ylabel('Coherence [-]');
            if exist('sgtitle') ~= 0, sgtitle(sprintf('Shaker point %d | Acc %d-%d', IP, I, I+2)); end
            drawnow;

            accepted = strcmpi(input('  Accept? (y/n): ', 's'), 'y');
            if ~accepted
                fprintf('  Record rejected — repeat.\n');
            end
        end

        % Store FRFs and coherences (columns = accelerometer positions)
        for k = 0:2
            Data.H{IP, I+k}   = H(:, k+1);
            Data.Coh{IP, I+k} = Co(:, k+1);
        end

        fname = sprintf('shaker_IP%02d_Acc%02d-%02d', IP, I, I+2);
        save(fullfile(SAVE_PATH, [fname '.mat']), 'acc', 'force', 'H', 'Co', 'freqh');
        saveas(gcf, fullfile(SAVE_PATH, [fname '.png']));
        fprintf('  Saved: %s\n', fname);
    end
end

Data.freqData = freqh;
save(fullfile(SAVE_PATH, 'shaker_all_points.mat'), 'Data');
fprintf('\nMeasurement complete. Data saved to: %s\n', SAVE_PATH);
