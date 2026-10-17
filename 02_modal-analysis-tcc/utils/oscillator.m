function wave = oscillator(wavetype, duration, frequency, gate, phase, sample_rate)
% =============================================================================
% oscillator.m — Test-signal generator (sine, sweeps, white/pink noise)
% =============================================================================
% THIRD-PARTY CODE (simplified)
%   Original author : W. Owen Brimijoin — MRC Institute of Hearing Research
%   Source          : MATLAB Central File Exchange — "oscillator" (v1.5, 2014)
%   Changes         : reduced to the waveforms used in this project
%                     (sinusoid, linear/log sweep, white and pink noise);
%                     input handling replaced by default arguments.
%
% Usage:
%   wave = oscillator(wavetype, duration, frequency, gate, phase, sample_rate)
%
% Inputs:
%   wavetype    : 'Sinusoid' | 'Linear Sweep' | 'Log Sweep' |
%                 'White Noise' | 'Pink Noise'          (case insensitive)
%   duration    : signal length [s]                     (default 1)
%   frequency   : [Hz]; [f_start f_stop] for sweeps; ignored for noise
%                                                        (default 440)
%   gate        : raised-cosine on/off ramp duration [s] (default 0)
%   phase       : starting phase [rad]                   (default 0)
%   sample_rate : [Hz]                                   (default 44100)
%
% Output:
%   wave        : column vector scaled to a peak absolute value of 1
%
% Example (shaker excitation, see shaker_measurement.m):
%   xs = oscillator('White Noise', 41, [], 0.1, 0, 5120);
% =============================================================================

if nargin < 1 || isempty(wavetype),    wavetype    = 'Sinusoid'; end
if nargin < 2 || isempty(duration),    duration    = 1;          end
if nargin < 3 || isempty(frequency),   frequency   = 440;        end
if nargin < 4 || isempty(gate),        gate        = 0;          end
if nargin < 5 || isempty(phase),       phase       = 0;          end
if nargin < 6 || isempty(sample_rate), sample_rate = 44100;      end

if 2*gate > duration
    error('oscillator: ''gate'' must be at most half the duration.');
end

num_samples = floor(sample_rate * duration);
t = linspace(0, duration, num_samples)';

switch lower(wavetype)
    case 'sinusoid'
        wave = sin(phase + 2*pi*frequency(1)*t);

    case {'linear sweep', 'log sweep'}
        if numel(frequency) ~= 2
            error('oscillator: sweeps require frequency = [f_start f_stop].');
        end
        if strcmpi(wavetype, 'linear sweep')
            f_inst = linspace(frequency(1), frequency(2), num_samples)';
        else
            f_inst = logspace(log10(frequency(1)), log10(frequency(2)), num_samples)';
        end
        % Integrate instantaneous frequency to obtain the phase
        wave = sin(phase + 2*pi*cumsum(f_inst)/sample_rate);

    case {'white noise', 'white'}
        wave = 2*(rand(num_samples, 1) - 0.5);

    case {'pink noise', 'pink'}
        % 1/f power spectrum shaped in the frequency domain
        f = [linspace(0, .5, floor(num_samples/2)), ...
             fliplr(linspace(0, .5, ceil(num_samples/2)))]';
        mag = 1 ./ sqrt(f);
        mag(isinf(mag)) = 0;
        wave = real(ifft(mag .* exp(2i*pi*rand(num_samples, 1))));

    otherwise
        error('oscillator: unrecognized waveform type ''%s''.', wavetype);
end

wave = wave ./ max(abs(wave));              % Normalize to +/- 1

if gate > 0                                 % Raised-cosine on/off ramps
    ramp = cos(linspace(0, pi/2, floor(sample_rate*gate))').^2;
    n = length(ramp);
    wave(1:n)         = wave(1:n)         .* flipud(ramp);
    wave(end-n+1:end) = wave(end-n+1:end) .* ramp;
end
end
