% =============================================================================
% modal_parameter_estimator.m — Frequency, Damping, and Mode Shape Extraction
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          Undergraduate Thesis (TCC) — UFSM
%
% Description:
%   Extracts modal parameters (natural frequencies, damping ratios, mode shapes)
%   from the measured FRF matrix using the stabilization diagram approach
%   via modalsd() (MATLAB Signal Processing Toolbox, R2017a+).
%
%   Processing pipeline:
%     Hvec [nlines x n_FRF] → CMIF → modalsd() stabilization → peak selection
%       → mode shapes (CMIF columns) → MAC matrix → animate_mode()
%
% Requirements:
%   - Signal Processing Toolbox (modalsd) — MATLAB only; in Octave the
%     stabilization diagram is skipped and fn_modes is set by hand
%     (use CMIF.m + pickpeaks.m to locate the resonances).
%   - Hvec and freqh loaded in workspace (from hammer or shaker measurement)
%
% Inputs (workspace variables):
%   Hvec  - FRF matrix [nlines x n_FRF]  (columns = individual FRFs)
%   freqh - Frequency vector [Hz]
%   Fs    - Sampling rate [Hz]
%
% Outputs (workspace variables):
%   fn_modes   - Natural frequencies of selected modes [1 x n_modes]
%   zeta_modes - Damping ratios of selected modes      [1 x n_modes]
%   Phi        - Mode shape matrix [n_DOF x n_modes]
% =============================================================================

%% =========================================================================
%  SECTION 1 — LOAD DATA
%  Load your own measurement file (not included in the repository)
% ==========================================================================

% load(fullfile('data', 'frf_measurement.mat'));   % your file: Hvec [nlines x 100], freqh, Fs

% Verify required variables
required = {'Hvec', 'freqh', 'Fs'};
for k = 1:length(required)
    if ~exist(required{k}, 'var')
        error('Variable "%s" not found. Load your data file first.', required{k});
    end
end

%% =========================================================================
%  SECTION 2 — FRF MATRIX SORTING
%  Separates FRFs by accelerometer group (X, Y, Z axes)
%  for a 10×10 measurement grid (100 impact-accelerometer combinations)
% ==========================================================================

nlines = size(Hvec, 1);

% Indices corresponding to axis columns (1,2,3 per impact point → X,Y,Z)
idx_xyz = false(1, 100);
for p = 1:100
    col = mod(p-1, 10) + 1;   % Column in 10×10 grid
    if col <= 3
        idx_xyz(p) = true;
    end
end

Hvec_xyz = Hvec(:, idx_xyz);    % 30 columns: 3 axes × 10 points

%% =========================================================================
%  SECTION 3 — STABILIZATION DIAGRAM (modalsd, MATLAB Signal Processing Toolbox)
% ==========================================================================

MAX_MODES = 35;     % Maximum number of modes to search

% Display stabilization diagram
has_modalsd = exist('modalsd') ~= 0;          % MATLAB Signal Processing Toolbox only
if has_modalsd
    figure('Name', 'Stabilization Diagram');
    fn_all = modalsd(Hvec_xyz, freqh, Fs, 'MaxModes', MAX_MODES);
else
    warning('modalsd not available (Octave): skipping the stabilization diagram.');
end

%% =========================================================================
%  SECTION 4 — SELECT MODES FROM STABILIZATION DIAGRAM
%  Inspect the stabilization diagram and choose the row (model order)
%  that best shows stable poles at the resonance frequencies of interest.
% ==========================================================================

if has_modalsd
    selected_row = 34;              % Edit this after inspecting the diagram
    tempFreqs = fn_all(selected_row, :);
    tempFreqs = tempFreqs(~isnan(tempFreqs) & tempFreqs > 0);

    fprintf('Selected %d candidate frequencies:\n', length(tempFreqs));
    fprintf('  %.2f Hz\n', tempFreqs);
end

%% =========================================================================
%  SECTION 5 — EXTRACT MODAL PARAMETERS
%  Pick resonance frequencies visually from CMIF or use pickpeaks()
% ==========================================================================

% Example: manually selected frequencies [Hz]
% Edit this vector with your identified resonances
fn_modes = [125.4, 245.8, 387.2, 512.6, 680.3];  % [Hz] — replace with your values

n_modes = length(fn_modes);
fprintf('\nExtracting %d modes...\n', n_modes);

% For each mode, find nearest CMIF peak
[~, idx_modes] = arrayfun(@(f) min(abs(freqh - f)), fn_modes);

% Mode shapes by quadrature peak picking: imag(H) at each resonance
% (valid for lightly damped structures; abs(H) would lose the sign)
Phi = zeros(size(Hvec, 2), n_modes);
for m = 1:n_modes
    Phi(:, m) = imag(Hvec(idx_modes(m), :))';   % quadrature picking keeps the sign
    Phi(:, m) = Phi(:, m) / max(abs(Phi(:, m)));    % Normalize
end

fprintf('Mode frequencies extracted. Review Phi matrix for mode shapes.\n');

%% =========================================================================
%  SECTION 6 — MAC MATRIX
% ==========================================================================

MAC_mat = zeros(n_modes, n_modes);
for i = 1:n_modes
    for j = 1:n_modes
        MAC_mat(i,j) = abs(Phi(:,i)' * Phi(:,j))^2 / ...
                       ((Phi(:,i)' * Phi(:,i)) * (Phi(:,j)' * Phi(:,j)));
    end
end

figure('Name', 'MAC Matrix');
imagesc(MAC_mat); colorbar; colormap('jet');
caxis([0 1]);
tick_lbl = arrayfun(@(f) sprintf('%.0f Hz', f), fn_modes, 'UniformOutput', false);
set(gca, 'XTick', 1:n_modes, 'YTick', 1:n_modes, 'XTickLabel', tick_lbl, 'YTickLabel', tick_lbl);
title('Modal Assurance Criterion (MAC)');
xlabel('Mode'); ylabel('Mode');

%% =========================================================================
%  SECTION 7 — MODE SHAPE ANIMATION (optional)
% ==========================================================================

% Reshape mode shape to 10×10 measurement grid and animate
% mode_idx = 1;   % Select mode to animate
% Psi = reshape(Phi(1:100, mode_idx), 10, 10);
% animate_mode(Psi, 2, fn_modes(mode_idx));
