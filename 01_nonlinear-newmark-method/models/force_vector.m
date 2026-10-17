% =============================================================================
% force_vector.m — External Force Vector Generator (1-DOF)
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Generates the external force vector F [1 x Np] for a single-DOF system
%   based on the selected excitation type (ccase).
%
% Inputs:
%   t      - Time vector                    [1 x Np]
%   ccase  - Excitation type (1 to 5):
%              1 = No force (free vibration)
%              2 = Unit impulse (duration = dt)
%              3 = Constant step force
%              4 = Harmonic (sinusoidal) force
%              5 = Step force that ends at t = Tf/2
%   Np     - Number of time steps           [integer]
%   exfreq - Excitation frequency           [Hz]  (used for ccase = 4)
%   Fmag   - Force magnitude                [N]
%   timp   - Impulse start time             [s]   (used for ccase = 2)
%   timpf  - Impulse end time               [s]   (used for ccase = 2)
%   dt     - Time step                      [s]
%   Ne     - Number of DOFs                 [integer]
%
% Outputs:
%   F      - Force matrix                   [Ne x Np]
%   force  - Scalar force amplitude used    [N]
%
% =============================================================================

function [F, force] = force_vector(t, ccase, Np, exfreq, Fmag, timp, timpf, dt, Ne)

F     = zeros(Ne, Np);
force = Fmag;

switch ccase

    case 1
        % ---- Free vibration (no external force) ----
        % F remains zero
        F = repmat(0, Ne, Np);

    case 2
        % ---- Unit impulse ----
        % Amplitude scaled so that integral(F)*dt = Fmag (unit impulse approximation)
        if timp == 0
            force = 2 * (Fmag / dt);
        else
            force = Fmag / dt;
        end
        for k = 1:Np
            if t(k) >= timp && t(k) <= timpf
                F(:,k) = force;
            end
        end

    case 3
        % ---- Constant step force ----
        F = repmat(Fmag, Ne, Np);

    case 4
        % ---- Harmonic force: F = Fmag * sin(2*pi*exfreq*t) ----
        F = Fmag * sin(2*pi * exfreq * t);

    case 5
        % ---- Step force ending at half-time ----
        halfNp = round(Np/2);
        F(:, 1:halfNp)       = Fmag;
        F(:, halfNp+1:end)   = 0;

    otherwise
        error('force_vector: unknown excitation case ccase = %d.', ccase);
end

end
