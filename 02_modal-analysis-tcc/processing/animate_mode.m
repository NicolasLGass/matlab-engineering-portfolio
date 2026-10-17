% =============================================================================
% animate_mode.m — 3D Mode Shape Animation
% =============================================================================
% THIRD-PARTY CODE (adapted)
%   Original author : unknown (snippet from experimental modal analysis
%                     teaching material)
%   Adapted by      : Nicolas Londero Gass — TCC UFSM
%   Changes         : VideoWriter export (aviobj/addframe removed), input
%                     defaults, normalization, MATLAB/Octave portability
%
% Description:
%   Animates a mode shape given as a 2D matrix (measurement grid).
%   The mode shape oscillates as Psi * sin(phase), cycling through a
%   full oscillation over the specified animation duration.
%
%   Optionally exports the animation to an .avi video file.
%
% Inputs:
%   Psi      - Mode shape matrix (reshaped to measurement grid) [rows x cols]
%   Maxtime  - Number of oscillation cycles to animate          [integer]
%   Freq     - Modal frequency for title display [Hz]           (optional)
%   Name     - Output filename (without extension) for .avi export
%              Leave empty ('') or omit to display only.
%
% Outputs:
%   y        - Always 0 (function used for side-effect animation only)
%
% Usage:
%   % Display only
%   animate_mode(Psi_mode1, 3, 125.4);
%
%   % Export to video
%   animate_mode(Psi_mode1, 3, 125.4, 'mode1_animation');
%
% =============================================================================

function y = animate_mode(Psi, Maxtime, Freq, Name)

%% --- Input validation ---
if nargin < 2
    Maxtime = 2;
end
if nargin < 3 || isempty(Freq)
    Freq = NaN;
    freq_str = 'Unknown';
else
    freq_str = sprintf('%.2f Hz', Freq);
end

export_video = false;
if nargin >= 4 && ~isempty(Name) && ~strcmp(Name, 'NaN')
    export_video = true;
end

%% --- Normalize mode shape ---
Psi = squeeze(Psi);
Psi = Psi ./ max(max(abs(Psi)));

%% --- Set up video writer (if exporting) ---
if export_video && exist('VideoWriter') == 0
    warning('VideoWriter not available (Octave without the video package): animation only.');
    export_video = false;
end
if export_video
    vid = VideoWriter([Name '.avi']);
    open(vid);
    fprintf('Exporting animation "%s.avi" at %s ...\n', Name, freq_str);
end

%% --- Animation loop ---
phase_vec = -pi : 0.1 : Maxtime * pi;
n_frames  = length(phase_vec);

figure(1); clf;

for k = 1:n_frames
    surf(Psi .* sin(phase_vec(k)));
    colormap('hot');
    axis equal; shading interp;
    zlim([-1 1]); caxis([-1 1]);
    title(sprintf('Mode Shape Animation: %s', freq_str), 'FontSize', 12);
    xlabel('Node X'); ylabel('Node Y'); zlabel('Normalized amplitude');
    drawnow;

    if export_video
        frame = getframe(gca);
        writeVideo(vid, frame);
    end
end

if export_video
    close(vid);
    fprintf('Video saved: %s.avi\n', Name);
end

y = 0;
end
