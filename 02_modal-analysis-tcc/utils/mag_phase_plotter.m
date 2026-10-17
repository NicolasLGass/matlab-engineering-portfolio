% =============================================================================
% mag_phase_plotter.m — FRF Magnitude and Phase Plots
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          Undergraduate Thesis (TCC) — UFSM
%
% Description:
%   Generates formatted magnitude and phase plots for a set of FRFs stored
%   in a cell array H{IP, acc}, where IP = impact point and acc = channel.
%
%   Two plot types are available:
%     1. Full 10×10 grid (all combinations)
%     2. Individual Bode-style plots for selected FRFs
%
% Inputs (workspace variables expected):
%   H     - FRF cell array {IP, acc} — complex [nlines x 1] per cell
%   freqh - Frequency vector [Hz]
%   XI    - Frequency axis lower limit [Hz]  (e.g., 30)
%   XF    - Frequency axis upper limit [Hz]  (e.g., 1000)
%   db    - Flag: 1 = plot in dB, 0 = plot in linear scale [m/s²/N]
%
% Usage:
%   Set the parameters below, load your H cell array, and run.
% =============================================================================

%% =========================================================================
%  PARAMETERS (edit here)
% ==========================================================================

XI = 30;        % Frequency axis lower limit [Hz]
XF = 1000;      % Frequency axis upper limit [Hz]
db = 1;         % 1 = dB, 0 = linear (m/s²/N)

n_IP  = 10;     % Number of impact points
n_acc = 10;     % Number of accelerometer channels

%% =========================================================================
%  PLOT 1 — Full 10×10 FRF Grid (Magnitude)
% ==========================================================================

figure('Name', 'FRF Magnitude Matrix (10x10)');
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

n = 1;
for IP = 1:n_IP
    for acc = 1:n_acc
        subplot(n_IP, n_acc, n);
        if db
            plot(freqh, mag2db(abs(H{IP, acc})), 'LineWidth', 0.8);
            ylabel('dB');
        else
            semilogy(freqh, abs(H{IP, acc}), 'LineWidth', 0.8);
            ylabel('[m/s²/N]');
        end
        xlim([XI XF]);
        title(sprintf('H(%d,%d)', IP, acc), 'FontSize', 7);
        set(gca, 'FontSize', 6);
        n = n + 1;
    end
end
if exist('sgtitle') ~= 0, sgtitle('FRF Magnitude Matrix', 'FontSize', 12); end

%% =========================================================================
%  PLOT 2 — Full 10×10 FRF Grid (Phase)
% ==========================================================================

figure('Name', 'FRF Phase Matrix (10x10)');
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

n = 1;
for IP = 1:n_IP
    for acc = 1:n_acc
        subplot(n_IP, n_acc, n);
        plot(freqh, angle(H{IP, acc}) * 180/pi, 'LineWidth', 0.8);
        xlim([XI XF]); ylim([-190 190]);
        ylabel('°'); title(sprintf('H(%d,%d)', IP, acc), 'FontSize', 7);
        set(gca, 'FontSize', 6);
        n = n + 1;
    end
end
if exist('sgtitle') ~= 0, sgtitle('FRF Phase Matrix', 'FontSize', 12); end

%% =========================================================================
%  PLOT 3 — Bode Diagram (select specific FRF)
% ==========================================================================

IP_sel  = 1;    % Select impact point
acc_sel = 1;    % Select accelerometer channel

figure('Name', sprintf('Bode — H(%d,%d)', IP_sel, acc_sel));

subplot(2,1,1);
if db
    plot(freqh, mag2db(abs(H{IP_sel, acc_sel})), 'b', 'LineWidth', 1.5);
    ylabel('Magnitude [dB]');
else
    semilogy(freqh, abs(H{IP_sel, acc_sel}), 'b', 'LineWidth', 1.5);
    ylabel('Magnitude [m/s²/N]');
end
grid on; xlim([XI XF]);
title(sprintf('FRF H_{%d,%d}', IP_sel, acc_sel));

subplot(2,1,2);
plot(freqh, angle(H{IP_sel, acc_sel}) * 180/pi, 'r', 'LineWidth', 1.5);
grid on; xlim([XI XF]); ylim([-190 190]);
xlabel('Frequency [Hz]'); ylabel('Phase [°]');
