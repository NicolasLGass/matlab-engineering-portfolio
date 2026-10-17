function mac = MAC(phi, show_plot)
% =============================================================================
% MAC.m — Modal Assurance Criterion matrix
% =============================================================================
% THIRD-PARTY CODE (simplified)
%   Original author : Mohammad Farshchin (University of Manitoba), 2015
%   Source          : MATLAB Central File Exchange — "MAC"
%   Changes         : loop replaced by vectorized expression, optional plot,
%                     complex mode shapes supported (conjugate transpose).
%
% Description:
%   MAC(i,j) = |phi_i' * phi_j|^2 / ((phi_i' * phi_i) * (phi_j' * phi_j))
%   Values close to 1 indicate correlated mode shapes; close to 0, distinct.
%   Reference: Allemang & Brown (1982), Proc. IMAC I, 110–116.
%
% Inputs:
%   phi       : [n_dof x n_modes] mode shape matrix (real or complex)
%   show_plot : (optional) plot the matrix as 3-D bars (default: true)
%
% Output:
%   mac       : [n_modes x n_modes] MAC matrix
%
% Example:
%   mac = MAC(Phi);
% =============================================================================

if nargin < 2
    show_plot = true;
end

cross = abs(phi' * phi).^2;                 % |phi_i' * phi_j|^2
auto  = real(diag(phi' * phi));             % phi_i' * phi_i
mac   = cross ./ (auto * auto.');

if show_plot
    figure;
    if exist('bar3') ~= 0
        bar3(mac); zlim([0 1]);              % MATLAB
    else
        imagesc(mac); colorbar; axis square; caxis([0 1]);   % Octave has no bar3
        set(gca, 'XTick', 1:size(mac,1), 'YTick', 1:size(mac,1));
    end
    title('Modal Assurance Criterion (MAC)');
    xlabel('Mode j'); ylabel('Mode i');
end
end
