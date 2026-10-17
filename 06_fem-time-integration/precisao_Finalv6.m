%% Elementos Finitos 3 - Análise de precisao dos métodos de integração
% Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
%
% Elongacao de periodo (PE) e decaimento de amplitude (AD) em funcao de dt/T
% para vibracao livre de 1 GDL nao amortecido, u(0) = 1 (Bathe, Fig. 9.3):
%   PE = (T_numerico - T)/T           [%]  (T_numerico pelos cruzamentos por zero)
%   AD = perda de amplitude por ciclo [%]  (razao entre picos sucessivos)
clear variables; close all; clc;
warning off;

%% Dados

Tf= 0.9999;                     % Periodo natural T [s]
omg= 2*pi/Tf;                   % Frequencia natural [rad/s]
dt= 0.02:0.001:0.2;             % Passos de tempo -> dt/T de 0.02 a 0.2
Ncyc= 20;                       % Numero de periodos simulados por dt

M= 1; K= omg^2; C= 0;           % Sistema de 1 GDL

u0= 1; u0p= 0; u0pp= -omg^2;    % Condicoes iniciais

tht= 1.4;                       % Theta (Wilson)
delt= 0.5;                      % Delta (Newmark, delt>=0.5)
alp= 0.25;                      % Alpha (Newmark, alp>=0.25*(0.5+delt)^2)

metodos= {'Central diff.', 'Houbolt', 'Wilson \theta=1.4', 'Newmark (\delta=0.5, \alpha=0.25)'};
nm= numel(metodos);
PE= nan(nm, numel(dt));
AD= nan(nm, numel(dt));

%% Integracao e medicao de periodo e amplitude

tic
for ii= 1:numel(dt)
    Np= ceil(Ncyc*Tf/dt(ii));
    R= zeros(1, Np+1);                                  % Vibracao livre
    U= zeros(nm, Np+1);
    U(1,:)= MDC2(M,K,C,R,u0,u0p,u0pp,Np,dt(ii));        % Diferenca Central
    U(2,:)= MH3(M,K,C,R,u0,u0p,u0pp,dt(ii),Np);         % Houbolt
    U(3,:)= MWT1(M,K,C,R,u0,u0p,u0pp,dt(ii),Np,tht);    % Wilson theta
    U(4,:)= MNK1(M,K,C,R,u0,u0p,u0pp,dt(ii),Np,delt,alp);   % Newmark
    t= (0:Np)*dt(ii);

    for m= 1:nm
        u= U(m,:);
        % Periodo: cruzamentos ascendentes por zero, com interpolacao linear
        k= find(u(1:end-1) < 0 & u(2:end) >= 0);
        if numel(k) >= 2
            tz= t(k) - u(k).*(t(k+1)-t(k))./(u(k+1)-u(k));
            PE(m,ii)= (mean(diff(tz)) - Tf)/Tf*100;
        end
        % Amplitude: picos positivos (maximos locais)
        p= find(u(2:end-1) > u(1:end-2) & u(2:end-1) >= u(3:end)) + 1;
        p= p(u(p) > 0);
        if numel(p) >= 2
            AD(m,ii)= (1 - (u(p(end))/u(p(1)))^(1/(numel(p)-1)))*100;
        end
    end
end
toc;

%% Resultados

fprintf('dt/T     ');  fprintf('%-14s', 'PE DC', 'PE Houb', 'PE Wil', 'PE Newm'); fprintf('\n');
for ii= [1 31 81 131 181]
    fprintf('%.3f   ', dt(ii)/Tf); fprintf('%+8.3f %%    ', PE(:,ii)); fprintf('\n');
end

cores= {'r', 'm', 'b', 'k'};

figure(1)
for m= 1:nm
    plot(dt/Tf, PE(m,:), cores{m}, 'linewidth', 1.5); hold on;
end
grid on; set(gca, 'fontsize', 12);
xlabel('\Deltat/T'); ylabel('Period elongation (T_{num}-T)/T [%]');
title('Period elongation of the integration methods');
legend(metodos, 'location', 'northwest'); xlim([0 0.2]);

figure(2)
for m= 1:nm
    plot(dt/Tf, AD(m,:), cores{m}, 'linewidth', 1.5); hold on;
end
grid on; set(gca, 'fontsize', 12);
xlabel('\Deltat/T'); ylabel('Amplitude decay per cycle [%]');
title('Amplitude decay of the integration methods');
legend(metodos, 'location', 'northwest'); xlim([0 0.2]);
