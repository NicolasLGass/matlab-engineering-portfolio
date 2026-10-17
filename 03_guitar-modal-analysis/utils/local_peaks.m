function [pks, locs] = local_peaks(x, min_dist, min_height, min_prom)
% =============================================================================
% local_peaks.m — Local maxima with minimum spacing, height and prominence
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Portable replacement for findpeaks(x, 'MinPeakDistance', d,
%   'MinPeakHeight', h, 'MinPeakProminence', p) that runs in both MATLAB and
%   Octave without the Signal Processing Toolbox / signal package.
%
%   - Peaks are local maxima (plateaus count once).
%   - Prominence: height of the peak above the higher of the two lowest
%     points between it and the nearest higher sample on each side (same
%     definition as findpeaks). Removes noise ripples riding on a crest.
%   - When two peaks are closer than min_dist samples, the higher one is
%     kept (same rule as findpeaks' MinPeakDistance).
%
% Inputs:
%   x          - Signal                                   [N x 1]
%   min_dist   - Minimum distance between peaks [samples] (default 1)
%   min_height - Minimum peak value                       (default -Inf)
%   min_prom   - Minimum peak prominence                  (default 0)
%
% Outputs:
%   pks  - Peak values (in time order)
%   locs - Peak indices (in time order)
% =============================================================================

if nargin < 2 || isempty(min_dist),   min_dist = 1;      end
if nargin < 3 || isempty(min_height), min_height = -Inf; end
if nargin < 4 || isempty(min_prom),   min_prom = 0;      end

x = x(:);
N = numel(x);
d = diff(x);
locs = find([false; d(1:end-1) > 0 & d(2:end) <= 0; false]);
locs = locs(x(locs) >= min_height);

% Prominence filter
if min_prom > 0 && ~isempty(locs)
    prom = zeros(size(locs));
    for k = 1:numel(locs)
        p = locs(k);
        il = p - 1;                                 % walk left to a higher sample
        while il >= 1 && x(il) <= x(p), il = il - 1; end
        ir = p + 1;                                 % walk right to a higher sample
        while ir <= N && x(ir) <= x(p), ir = ir + 1; end
        left_min  = min(x(max(il, 1):p));
        right_min = min(x(p:min(ir, N)));
        prom(k) = x(p) - max(left_min, right_min);
    end
    locs = locs(prom >= min_prom);
end

% Enforce minimum distance: keep the highest peaks first
if min_dist > 1 && numel(locs) > 1
    [~, order] = sort(x(locs), 'descend');
    keep = false(size(locs));
    taken = false(size(x));
    for k = order(:)'
        p = locs(k);
        if ~taken(p)
            keep(k) = true;
            taken(max(1, p-min_dist+1) : min(N, p+min_dist-1)) = true;
        end
    end
    locs = locs(keep);
end

pks = x(locs);
end
