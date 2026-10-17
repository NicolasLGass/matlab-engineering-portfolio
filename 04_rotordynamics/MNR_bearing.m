%% Metodo de Newton-Raphson - equilibrio do mancal hidrodinamico curto
% Resolve o equilibrio de forcas do mancal para a excentricidade E e o
% angulo de atitude phi, com Jacobiano analitico (sem toolbox simbolico):
%   f1 = Ft*sin(phi) + Fr*cos(phi) - F0 = 0
%   f2 = Ft*cos(phi) - Fr*sin(phi)      = 0
%   Ft = K0*pi*E/(4*(1-E^2)^(3/2)),  Fr = K0*E^2/(1-E^2)^2
% Variaveis de entrada:
% K0  - eta*omega*R*L^3/c^2 (fator de forca do mancal curto)
% F0  - Carga estatica no mancal [N]
% est - Valores iniciais [E; phi]
% Saidas:
% E, phi - Excentricidade relativa e angulo de atitude [rad]

function [E, phi] = MNR_bearing(K0, F0, est)

x= est(:);
for i= 1:31
    s= x - bearing_jacobian(x, K0) \ bearing_residual(x, K0, F0);
    err= abs(s - x);
    x= s;
    if max(err) <= 1e-7                     % todas as variaveis convergiram
        break
    end
end
if max(err) > 1e-7
    disp('Resultado não converge para os valores iniciais estabelecidos.')
end
E= x(1);
phi= x(2);
end

function r = bearing_residual(x, K0, F0)
E= x(1); ph= x(2);
Ft= K0*pi*E/(4*(1-E^2)^1.5);                % Forca tangencial
Fr= K0*E^2/(1-E^2)^2;                       % Forca radial
r= [Ft*sin(ph) + Fr*cos(ph) - F0;
    Ft*cos(ph) - Fr*sin(ph)];
end

function J = bearing_jacobian(x, K0)
E= x(1); ph= x(2);
Ft= K0*pi*E/(4*(1-E^2)^1.5);
Fr= K0*E^2/(1-E^2)^2;
dFt= K0*(pi/4)*(1+2*E^2)/(1-E^2)^2.5;       % dFt/dE
dFr= K0*2*E*(1+E^2)/(1-E^2)^3;              % dFr/dE
J= [dFt*sin(ph) + dFr*cos(ph),   Ft*cos(ph) - Fr*sin(ph);
    dFt*cos(ph) - dFr*sin(ph),  -Ft*sin(ph) - Fr*cos(ph)];
end
