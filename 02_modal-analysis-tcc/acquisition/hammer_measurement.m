% =============================================================================
% hammer_measurement.m — Impact Hammer Modal Measurement (NI cDAQ)
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          Undergraduate Thesis (TCC) — UFSM
%
% Description:
%   Real-time data acquisition routine for experimental modal analysis
%   using an impact hammer and three IEPE accelerometers connected to a
%   National Instruments cDAQ-9234 module.
%
%   Acquisition strategy:
%     - 10 impact points × 10 accelerometer groups (3 channels each)
%     - 6 averages per measurement point
%     - H1 FRF estimator computed per impact
%     - Bad impacts (double hits, low coherence) rejected interactively
%
%   Signal chain:
%     NI cDAQ-9234 (IEPE) → MATLAB DAQ Toolbox → FRF (H1) → Saved .mat
%
% Hardware:
%   - NI cDAQ-9234, 4 channels, 24-bit
%   - 3 × IEPE accelerometers (PCB / similar)
%   - 1 × force transducer (impact hammer)
%
% Requirements:
%   - MATLAB Data Acquisition Toolbox
%   - NI-DAQmx driver
%
% Usage:
%   1. Connect hardware and verify channel numbering.
%   2. Set SAVE_PATH to your target directory.
%   3. Run the script — follow on-screen prompts per measurement point.
% =============================================================================

clear; format long; clc; close all;
if exist('OCTAVE_VERSION', 'builtin')
    error('This acquisition script needs MATLAB + Data Acquisition Toolbox (NI-DAQmx); it does not run in Octave.');
end

%% =========================================================================
%  SECTION 1 — USER SETTINGS (edit here)
% ==========================================================================

SAVE_PATH = fullfile(pwd, 'data', 'hammer');    % Output directory for .mat files
Fs        = 5120;                               % Sampling rate [Hz]
%   Available rates for NI-9234: 12800, 6400, 5120, 4266, 3840, ...
%   Formula: (131072/256) / n, n = 1..31

T_rec    = 10;                                  % Record duration per impact [s]
N_avg    = 6;                                   % Number of valid averages to collect

NFFT = 1 * Fs;                                  % FFT length (= 1 s window)

%% =========================================================================
%  SECTION 2 — TRANSDUCER SENSITIVITIES
% ==========================================================================

sens_acc1  = 0.9964e-3;     % [V/(m/s²)] — Accelerometer channel 0 (X)
sens_acc2  = 1.036e-3;      % [V/(m/s²)] — Accelerometer channel 1 (Y)
sens_acc3  = 1.025e-3;      % [V/(m/s²)] — Accelerometer channel 2 (Z)
sens_force = 20.35e-3;      % [V/N]      — Force transducer channel 3

cal1 = 1 / sens_acc1;
cal2 = 1 / sens_acc2;
cal3 = 1 / sens_acc3;
cal4 = 1 / sens_force;

%% =========================================================================
%  SECTION 3 — DAQ SESSION INITIALIZATION
% ==========================================================================

if ~exist(SAVE_PATH, 'dir'), mkdir(SAVE_PATH); end    % Create output folder if needed

d = daq.getDevices;
s = daq.createSession('ni');
s.Rate             = Fs;
s.DurationInSeconds = T_rec;

% Add 4 IEPE channels (0=AccX, 1=AccY, 2=AccZ, 3=Force)
[ch, ~] = addAnalogInputChannel(s, 'cDAQ1Mod1', [0, 1, 2, 3], 'IEPE');

fprintf('DAQ initialized: Fs = %d Hz, T_rec = %d s, NFFT = %d\n', Fs, T_rec, NFFT);

%% =========================================================================
%  SECTION 4 — DATA STRUCTURE INITIALIZATION
% ==========================================================================

Data = struct();                                 % results container (H, Coh, settings)

Data.SampRate     = Fs;
Data.MeasuringTime = T_rec;
Data.Nbits        = 24;
Data.Sensibility.Acc1  = sens_acc1;
Data.Sensibility.Acc2  = sens_acc2;
Data.Sensibility.Acc3  = sens_acc3;
Data.Sensibility.Force = sens_force;

%% =========================================================================
%  SECTION 5 — MEASUREMENT LOOP
%  Outer loop: accelerometer groups (3 channels at a time, 10 groups total)
%  Inner loop: impact points (10 points per structure)
% ==========================================================================

for I = 1:3:10              % Accelerometer group index (channels 1-3, 4-6, ...)

    for IP = 1:10           % Impact point index

        TIP  = 10;          % Total impact points
        nn   = 1;           % Loop control flag
        n_ok = 0;           % Valid average counter

        fprintf('\n--- Group %d | Impact point %d/%d ---\n', I, IP, TIP);

        Hx_acc = zeros(NFFT/2+1, N_avg);
        Hy_acc = zeros(NFFT/2+1, N_avg);
        Hz_acc = zeros(NFFT/2+1, N_avg);
        Cx_acc = zeros(NFFT/2+1, N_avg);
        Cy_acc = zeros(NFFT/2+1, N_avg);
        Cz_acc = zeros(NFFT/2+1, N_avg);

        while nn && n_ok < N_avg

            % Prompt user to impact
            fprintf('  Ready for impact %d/%d — press ENTER to arm\n', n_ok+1, N_avg);
            pause;

            % Acquire data
            raw = s.startForeground();

            % Apply calibration
            accx  = raw(:,1) * cal1;
            accy  = raw(:,2) * cal2;
            accz  = raw(:,3) * cal3;
            force = raw(:,4) * cal4;

            t = (0:length(accx)-1)' / Fs;

            % FRF estimation (H1) and coherence
            ov  = 0;            % No overlap for impact testing
            win = rectwin(NFFT);

            [Hx, ~]     = tfestimate(force, accx,  win, ov, NFFT, Fs);
            [Hy, ~]     = tfestimate(force, accy,  win, ov, NFFT, Fs);
            [Hz, freqh] = tfestimate(force, accz,  win, ov, NFFT, Fs);

            [Cx, ~]     = mscohere(force, accx,  win, ov, NFFT, Fs);
            [Cy, ~]     = mscohere(force, accy,  win, ov, NFFT, Fs);
            [Cz, ~]     = mscohere(force, accz,  win, ov, NFFT, Fs);

            % Quick-check plot for the user
            figure(99); clf;
            subplot(2,1,1); plot(t, force);
            xlabel('t [s]'); ylabel('Force [N]');
            title(sprintf('Impact check — Point %d | Avg %d/%d', IP, n_ok+1, N_avg));
            grid on;
            subplot(2,1,2); semilogy(freqh, abs(Hz));
            xlabel('f [Hz]'); ylabel('|H_z| [m/s²/N]'); grid on;
            drawnow;

            % Accept or reject
            resp = input('  Accept? (y/n): ', 's');
            if strcmpi(resp, 'y')
                n_ok = n_ok + 1;
                Hx_acc(:, n_ok) = Hx;
                Hy_acc(:, n_ok) = Hy;
                Hz_acc(:, n_ok) = Hz;
                Cx_acc(:, n_ok) = Cx;
                Cy_acc(:, n_ok) = Cy;
                Cz_acc(:, n_ok) = Cz;
            else
                fprintf('  Impact rejected — repeat.\n');
            end
        end

        % Average FRFs across valid impacts
        Data.H{IP, I}   = mean(Hx_acc, 2);
        Data.H{IP, I+1} = mean(Hy_acc, 2);
        Data.H{IP, I+2} = mean(Hz_acc, 2);

        Data.Coh{IP, I}   = mean(Cx_acc, 2);
        Data.Coh{IP, I+1} = mean(Cy_acc, 2);
        Data.Coh{IP, I+2} = mean(Cz_acc, 2);

        % Save after each impact point
        fname = fullfile(SAVE_PATH, sprintf('hammer_IP%02d_Gr%02d.mat', IP, I));
        freqh_out = freqh;  %#ok<NASGU>
        save(fname, 'Data', 'freqh_out');
        fprintf('  Saved: %s\n', fname);

    end % IP loop
end % I loop

fprintf('\nMeasurement complete. Data saved to: %s\n', SAVE_PATH);
