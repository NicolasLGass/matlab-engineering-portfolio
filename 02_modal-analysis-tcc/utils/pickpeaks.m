function [peaks, criterion] = pickpeaks(V, select, display)
% =============================================================================
% pickpeaks.m — Scale-space peak picking
% =============================================================================
% THIRD-PARTY CODE (simplified)
%   Original author : Antoine Liutkus — Inria Nancy, 2015 (BSD-3-Clause)
%   Reference       : A. Liutkus, "Scale-Space Peak Picking", Research Report,
%                     Inria, 2015. https://hal.inria.fr/hal-01103123
%   Changes         : comments condensed, unused variables removed,
%                     histc() kept (runs in both MATLAB and Octave).
%                     Algorithm unchanged.
%
% Used here to pick resonance candidates from |FRF| or CMIF curves.
%
% Inputs:
%   V       : data vector
%   select  : > 1      -> number of peaks to return
%             0 < x < 1 -> criterion threshold (closer to 1 = fewer peaks)
%   display : (optional) plot data, criterion and peaks (default 0)
%
% Outputs:
%   peaks     : indices of the detected peaks
%   criterion : per-sample peak likelihood, normalized to [0, 1]
%
% Algorithm: the data are smoothed with increasing horizons; local maxima
% that survive across scales accumulate a score weighted by the scale.
% Points with the largest accumulated score are returned as peaks.
% ---------------------------------------------------------------------
% Copyright (C) 2015, Inria, Antoine Liutkus
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met:
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Inria nor the names of its contributors may
%       be used to endorse or promote products derived from this software
%       without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL INRIA BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% =============================================================================

if nargin < 3, display = 0; end
if nargin < 2, select  = 0; end

V = V(:) - min(V(:));
Vorig = V;
n = length(V);

buffer    = zeros(n,1);
criterion = zeros(n,1);
if select < 1
    minDist = n/20;
else
    minDist = n/select;
end
horizons = unique(round(logspace(0,2,50)/100*ceil(n/20)));

% Neighbour matrix used to find local maxima (faster than findpeaks)
tempMat = zeros(n,3);
tempMat(1,1)   = inf;
tempMat(end,3) = inf;

for is = 1:length(horizons)

    % Smooth with a half-sine kernel (FFT-based convolution)
    horizon = horizons(is);
    if horizon > 1
        w = max(eps, sin(2*pi*(0:(horizon-1))/2/(horizon-1)));
        w = w/sum(w);
        V = real(ifft(fft(V(:),n+horizon).*fft(w(:),n+horizon)));
        V = V(1+floor(horizon/2):end-ceil(horizon/2));
    end

    % Local maxima of the smoothed data
    tempMat(2:end,1)   = V(1:end-1);
    tempMat(:,2)       = V(:);
    tempMat(1:end-1,3) = V(2:end);
    [~, posMax] = max(tempMat, [], 2);
    I = find(posMax == 2)';

    newBuffer = zeros(size(buffer));
    if is == 1
        newBuffer(I) = Vorig(I);            % keep all maxima at the finest scale
    else
        old = find(buffer)';
        if isempty(old)
            continue;
        end
        % Link each new maximum to the closest maximum of the previous scale
        % (nearest-neighbour trick by Roger Stafford, MATLAB Central)
        [c, p] = sort(old);
        [~, ic] = histc(I, [-inf, (c(1:end-1)+c(2:end))/2, inf]);   % histc: MATLAB and Octave
        iOld = p(ic);
        d    = abs(I - old(iOld));
        neighbours = iOld(d < minDist);
        if ~isempty(neighbours)
            newBuffer(old(neighbours)) = V(old(neighbours))*is^2;
        end
    end
    buffer    = newBuffer;
    criterion = criterion + newBuffer;
end

criterion = criterion/max(criterion);

if select < 1
    peaks = find(criterion > select);
else
    [~, order] = sort(criterion, 'descend');
    peaks = order(1:select);
end

if display
    clf;
    plot(Vorig, 'LineWidth', 2); hold on;
    plot(criterion*max(Vorig), 'r');
    plot(peaks, Vorig(peaks), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    grid on;
    title('Scale-space peak detection');
    legend('data', 'computed criterion', 'selected peaks');
end
end
