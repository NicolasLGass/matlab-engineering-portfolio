% =============================================================================
% find_decay.m — Free-Decay Envelope Extraction
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%          MEAV — UFSM (2019)
%
% Description:
%   Extracts the upper decay envelope from a vibration free-decay signal.
%   Steps:
%     1. Remove signal delay (signal starts at its absolute maximum)
%     2. Extract upper envelope using peak detection with spacing n
%
%   Used as the first step in logarithmic decrement damping estimation.
%
% Inputs:
%   x    - Input vibration signal (time-domain)             [N x 1]
%   n    - Envelope smoothing spacing (samples between peaks) [integer]
%            Larger n = smoother envelope, fewer points
%            Typical: 50–400 depending on signal frequency and Fs
%   c    - (Optional) Crop signal to first c samples after delay onset
%   ab   - (Optional) 1 = force first peak positive (flip sign if needed)
%
% Outputs:
%   upper - Upper decay envelope (peak values)              [variable length]
%   fsig  - Signal after delay removal (same length as envelope input)
%
% Example:
%   [env, sig] = find_decay(accz, 180);
%   t_env = linspace(0, length(env)/Fs, length(env));
%   plot(t_env, env);
%
% =============================================================================

function [upper, fsig] = find_decay(x, n, c, ab)

narginchk(1, 4);
nargoutchk(0, 2);

%% --- Remove signal delay: start the decay at the absolute maximum ---
x = x(:);
[~, i0] = max(abs(x));

if nargin > 2 && ~isempty(c)
    dx = x(i0 : min(i0 + c - 1, end));
else
    dx = x(i0 : end);
end

%% --- Optional: force positive initial peak ---
if nargin > 3 && ab == 1
    if abs(min(dx)) > max(dx)
        dx = -dx;
    end
end

fsig = dx;

%% --- Extract upper envelope ---
% Spline through peaks at least n samples apart (equivalent to
% envelope(dx, n, 'peak'), portable to Octave)
[~, locs] = local_peaks(dx, n);
knots = unique([1; locs; numel(dx)]);
upper = interp1(knots, dx(knots), (1:numel(dx))', 'spline');

end
