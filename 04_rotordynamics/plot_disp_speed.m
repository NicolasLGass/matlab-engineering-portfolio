function plot_disp_speed(t, u, rpm, ylims, ylab, ttl)
% Grafico com dois eixos y: deslocamento (esquerda) e rotacao (direita).
% Usa plotyy, disponivel no MATLAB e no Octave (yyaxis so existe no MATLAB).
% t     - vetor de tempo [s]
% u     - deslocamento [m]
% rpm   - rotacao [rpm]
% ylims - limites do eixo de deslocamento [m]
% ylab  - rotulo do eixo de deslocamento
% ttl   - (opcional) titulo do subplot

[ax, h1, h2] = plotyy(t, u, t, rpm);
set(h1, 'Color', 'k', 'LineWidth', 0.85);
set(h2, 'Color', 'b', 'LineWidth', 1);
ylim(ax(1), ylims);
set(ax(1), 'FontSize', 10, 'XColor', 'k', 'YColor', 'k', ...
    'YTickMode', 'auto', 'XGrid', 'on', 'YGrid', 'on');
set(ax(2), 'FontSize', 10, 'XColor', 'k', 'YColor', 'b', 'YTickMode', 'auto');
xlabel(ax(1), 'Time [s]', 'FontWeight', 'bold');
ylabel(ax(1), ylab, 'FontWeight', 'bold');
ylabel(ax(2), 'Shaft speed [rpm]', 'FontWeight', 'bold');
if nargin > 5
    title(ax(1), ttl, 'FontSize', 12);
end
end
