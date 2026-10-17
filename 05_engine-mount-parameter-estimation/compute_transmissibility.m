function [Tf_all, fr_all] = compute_transmissibility(D, L, E, eta, M, w)
% =============================================================================
% compute_transmissibility.m — Rubber mount transmissibility and resonance
% =============================================================================
% Author : Nicolas Londero Gass
%          nicolas.gass@eac.ufsm.br
%
% Computes the force transmissibility and resonance frequency of a mass M
% on 4 identical cylindrical rubber mounts (complex stiffness, loss factor eta).
% Kept in its own file (instead of a local function) so that
% engine_mount_analysis.m also runs in GNU Octave.
%
% Inputs:
%   D, L, E - Mount diameter [m], length [m] and Young's modulus [Pa] (vectors)
%   eta     - Loss factor [-]
%   M       - Supported mass [kg]
%   w       - Angular frequency vector [rad/s]
% Outputs:
%   Tf_all  - |T(w)| for each configuration   [n x length(w)]
%   fr_all  - Resonance frequency [Hz]          [1 x n]
% =============================================================================

    n = length(D);
    Tf_all = zeros(n, length(w));
    fr_all = zeros(1, n);

    for i = 1:n
        A  = pi * D(i)^2 / 4;                      % Contact area [m²]
        Ki = (E(i) * (1 + 1i*eta) * A) / L(i);     % Single mount stiffness
        Kt = 4 * Ki;                                 % Total (4 mounts in parallel)

        Tf_all(i,:) = abs(Kt ./ (Kt - M * w.^2));
        fr_all(i)   = (1/(2*pi)) * sqrt(real(Kt)/M);
    end
end
