% =============================================================================
% validate_spring_pendulum.m — Nonlinear (elastic) pendulum, Géradin & Rixen
% =============================================================================
% Author : Nicolas Londero Gass
%          MSc Mechanical Engineering — UNICAMP (2025)
%          nicolas.gass@eac.ufsm.br
%
% Description:
%   Elastic pendulum example from Géradin & Rixen (2014), "Mechanical
%   Vibrations", Ch. 7, solved with the nonlinear Newmark predictor-corrector
%   scheme and Newton-Raphson iterations (MSc dissertation, appendix).
%
%   u = [u1; u2]: u1 vertical (upwards), u2 horizontal; pivot at the origin.
%       L = sqrt(u1^2 + u2^2)
%       M*u'' + k*(L - L0)/L * u = [-m*g; 0]
%
%   m = 1 kg, k = 30 N/m, g = 10 m/s^2, L0 = 1 m,
%   u(0) = [0; 1.5] m, u'(0) = 0,
%   gamma = 1/2, beta = 1/4, dt = 0.03 s, tolerance m*g*1e-5.
%
%   The book gives only plots, so the Newmark solution is also compared with
%   a high-accuracy ode45 reference, repeated for smaller time steps
%   (convergence order) and checked for total-energy conservation.
% =============================================================================

clear; close all; clc;

plotlentotrajetoria = 0;          % 1 = animate the trajectory
passoPlt = 20;                    % animation step

%% Parameters
M = [1 0;
     0 1];

g = 10;
L0 = 1;
err = 1*g*10^(-5);                % Tolerance

C = [0 0; 0 0];

L = @(u,v) sqrt((u^2)+(v^2));

k = 30*1;
K = @(u,v) [k*(L(u,v)-L0)*1/L(u,v)  0;
            0  k*(L(u,v)-L0)*1/L(u,v)];

% Tangent stiffness dP/du
Kt = @(u,v) [k*(1- L0/L(u,v)) + k*(L0*u^2)/(L(u,v))^3   k*((L0*u*v)/(L(u,v)^3));
             k*((L0*u*v)/(L(u,v)^3))   k*(1- L0/L(u,v)) + k*(L0*v^2)/(L(u,v))^3];

u0 = 0; v0 = 1.5;
u0p = 0; v0p = 0;

gamma = .5; beta = 1/4;

tf = 20;                          % Integration time

% Reference solution (ode45, tight tolerances)
odefun = @(tt, z) [z(3:4); M\([-M(1,1)*g; 0] - K(z(1),z(2))*z(1:2))];
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

dts = [0.03 0.01 0.003 0.001];    % 0.03 s = book time step
err5 = zeros(size(dts));
fprintf('  dt [s]   max error 0-5 s [m]   error > 0.1 m after   energy drift   Newton it/step\n');
fprintf('  -----------------------------------------------------------------------------------\n');

for idt = 1:numel(dts)
    dt = dts(idt);
    t = (0:dt:tf);                % Time vector
    Np = length(t);               % Number of time points

    val1 = -M*g;
    R = repmat(val1(:,1), 1, Np); % Gravity load

    %% Initial values
    Ne = size(R,1);               % Number of DOFs

    P0 = K(u0,v0)*[u0; v0];

    u0pp = M\(R(:,1)-P0);         % Initial acceleration

    a1 = (1-gamma)*dt;
    a2 = (.5-beta)*dt^2;
    a3 = (gamma/(beta*dt));
    a4 = (1/((dt^2)*beta));

    %% Newmark time integration
    n = 1;
    u = zeros(Ne, Np); up = zeros(Ne, Np); upp = zeros(Ne, Np);
    ue = zeros(Ne, Np+n); upe = zeros(Ne, Np+n); uppe = zeros(Ne, Np+n);
    Pe = zeros(Ne, Np+n); Psi = zeros(Ne, Np+n); dtu = zeros(Ne, Np+n);
    nit = zeros(1, Np-1);

    Ct = C;

    u(:,1) = [u0; v0];
    up(:,1) = [u0p; v0p];
    upp(:,1) = u0pp;

    for i = 1:Np-1
        % Predictor
        uppe(:,i+1) = 0;                                    % u..* - predicted acceleration
        upe(:,i+1) = up(:,i) + a1*upp(:,i);                 % u.*  - predicted velocity
        ue(:,i+1) = u(:,i) + dt*up(:,i) + a2*upp(:,i);      % u*   - predicted displacement

        contwhile = 0;
        convtest = 10;
        while norm(convtest) > err
            contwhile = contwhile+1;
            if contwhile > 100
                fprintf('No convergence after 100 iterations, i = %d\n', i);
                break;
            end

            uv = [ue(1,i+1); ue(2,i+1)];

            Pe(:,i+1) = K(ue(1,i+1),ue(2,i+1))*uv;
            Psi(:,i+1) = M*uppe(:,i+1) + Pe(:,i+1) - R(:,i+1);   % residual
            PsiOld = Psi(:,i+1);

            % Iteration matrix
            S = Kt(ue(1,i+1),ue(2,i+1)) + a3*Ct + a4*M;         % dPsi/du
            dtu(:,i+1) = -S\(Psi(:,i+1));                       % delta_u

            % Corrector
            u(:,i+1) = ue(:,i+1) + dtu(:,i+1);
            up(:,i+1) = upe(:,i+1) + a3*dtu(:,i+1);
            upp(:,i+1) = uppe(:,i+1) + a4*dtu(:,i+1);

            Pe(:,i+1) = K(ue(1,i+1),ue(2,i+1))*uv;
            Psi(:,i+1) = M*upp(:,i+1) + Pe(:,i+1) - R(:,i+1);   % residual

            convtest = Psi(:,i+1)-PsiOld;
            if norm(convtest) < err
                break
            end

            ue(:,i+1) = u(:,i+1);
            upe(:,i+1) = up(:,i+1);
            uppe(:,i+1) = upp(:,i+1);
        end
        nit(i) = contwhile;
    end

    %% Comparison with the reference solution
    [~, z] = ode45(odefun, t, [u0; v0; u0p; v0p], opts);
    e = max(abs(u - z(:,1:2).'), [], 1);
    err5(idt) = max(e(t <= 5));
    i01 = find(e > 0.1, 1);
    if isempty(i01), tdiv = '      -  '; else, tdiv = sprintf('%7.2f s', t(i01)); end
    E = 0.5*sum((M*up).*up, 1) + 0.5*k*(sqrt(sum(u.^2, 1)) - L0).^2 + M(1,1)*g*u(1,:);
    fprintf('  %6.3f        %9.2e             %s        %6.3f %%          %.2f\n', ...
            dt, err5(idt), tdiv, max(abs(E - E(1)))/abs(E(1))*100, mean(nit));

    if idt == 1                    % keep the book case for the plots
        tB = t; uB = u; upB = up; zB = z;
    end
end
pc = polyfit(log(dts), log(err5), 1);
fprintf('\n  Observed convergence order (0-5 s): %.2f  (theory: 2)\n', pc(1));

t = tB; u = uB; up = upB; z = zB;

%% Plots (book case, dt = 0.03 s; x = -u1 downwards, y = u2 as in the book)

figure(1);
plot(t,-u(1,:)); hold on; grid on;
plot(t,u(2,:), 'k--');
xlabel('Time [s]');
ylabel('Displacement [m]');
legend('x (vertical)','y (horizontal)');
title('Nonlinear pendulum - displacement (\Deltat = 0.03 s)');
ylim([-2 3])

figure(2);
plot(t,-up(1,:)); hold on; grid on;
plot(t,up(2,:), 'k--');
xlabel('Time [s]');
ylabel('Velocity [m/s]');
legend('dx/dt','dy/dt');
title('Nonlinear pendulum - velocity (\Deltat = 0.03 s)');

figure(3);
plot(u(2,:),u(1,:)); hold on; grid on;
plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
axis equal;
xlabel('y [m]');
ylabel('-x [m]');
legend('Trajectory', 'Pivot');
title('Nonlinear pendulum - trajectory');

figure(4);
plot(-u(1,:),-up(1,:)); hold on; grid on;
xlabel('x [m]');
ylabel('dx/dt [m/s]');
title('Nonlinear pendulum - phase plane in x');

figure(5);
plot(u(2,:),up(2,:)); hold on; grid on;
xlabel('y [m]');
ylabel('dy/dt [m/s]');
title('Nonlinear pendulum - phase plane in y');

figure(6);
plot(t, -z(:,1), 'k', 'LineWidth', 1.6); hold on; grid on;
plot(t, -u(1,:), '--', 'LineWidth', 1.2);
plot(t, z(:,2), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.6);
plot(t, u(2,:), ':', 'LineWidth', 1.4);
xlabel('Time [s]');
ylabel('Displacement [m]');
legend('x - ode45', 'x - Newmark', 'y - ode45', 'y - Newmark', 'Location', 'northeast');
title('Newmark (\Deltat = 0.03 s) vs ode45 reference');
xlim([0 10]); ylim([-2 3])

figure(7);
loglog(dts, err5, 'o-', 'LineWidth', 1.5); hold on; grid on;
loglog(dts, err5(1)*(dts/dts(1)).^2, 'k--');
xlabel('\Deltat [s]');
ylabel('Max. displacement error, 0-5 s [m]');
legend('Newmark vs ode45', 'slope 2', 'Location', 'northwest');
title('Time-step convergence');

%% Animated trajectory

if plotlentotrajetoria == 1
    figure(8);
    h = plot(u(2,1),u(1,1),'k'); hold('all'); grid on;
    xlabel('y [m]');
    ylabel('-x [m]');
    title('Nonlinear pendulum - trajectory');
    xlim([-2 2])
    ylim([-2.5 0.2])

    for i = (passoPlt+1):passoPlt:length(u(1,:))
        x = u(2,i-passoPlt:i);
        y = u(1,i-passoPlt:i);
        plot(x,y,'k');
        drawnow;
    end
end
