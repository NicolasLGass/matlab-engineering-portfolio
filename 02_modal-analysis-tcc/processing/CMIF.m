% =============================================================================
% CMIF.m — Complex Mode Indicator Function
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          Undergraduate Thesis (TCC) — UFSM
%
% Description:
%   Computes the Complex Mode Indicator Function (CMIF) from the FRF matrix.
%   The CMIF is the matrix of singular values of the FRF matrix [H(f)] at
%   each frequency line, computed via SVD:
%
%       [U, S, V] = svd(H(f_i))    for each frequency line i
%       CMIF(f_i) = diag(S)
%
%   Peaks in the CMIF curves indicate resonance frequencies.
%   The number of significant singular values at a peak indicates the
%   number of closely spaced or repeated modes.
%
% Inputs:
%   Hvec  - FRF matrix [nlines x n_acc x n_force]
%             nlines  = number of frequency lines
%             n_acc   = number of response (accelerometer) DOFs
%             n_force = number of excitation (force) DOFs
%   freqh - Frequency vector corresponding to rows of Hvec [nlines x 1]
%
% Outputs:
%   h     - CMIF matrix [nlines x n_force]  (singular values per frequency)
%   fig   - Handle to the generated plot
%
% Usage:
%   [h, fig] = CMIF(Hvec, freqh);
%
% Reference:
%   Allemang, R.J. & Brown, D.L. (1982). A Correlation Coefficient for
%   Modal Vector Analysis. Proc. IMAC I.
% =============================================================================

function [h, fig] = CMIF(Hvec, freqh)

nlines  = size(Hvec, 1);
n_acc   = size(Hvec, 2);
n_force = size(Hvec, 3);

h = zeros(nlines, n_force);

for i = 1:nlines
    % Singular values of the FRF slice at frequency line i
    sv = svd(squeeze(Hvec(i, :, :)));          % complex FRF matrix (keeps phase)
    h(i, :) = sv(1:n_force)';
end

% Plot
fig = figure('Name', 'CMIF — Complex Mode Indicator Function');
semilogy(freqh, h, 'LineWidth', 1.5);
grid on;
xlabel('Frequency [Hz]', 'FontSize', 11);
ylabel('Singular values [-]', 'FontSize', 11);
title('Complex Mode Indicator Function (CMIF)', 'FontSize', 12);
if n_force > 1
    leg = arrayfun(@(k) sprintf('SV %d', k), 1:n_force, 'UniformOutput', false);
    legend(leg, 'Location', 'northeast');
end

end
