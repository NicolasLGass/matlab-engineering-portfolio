%% Elementos Finitos 3 - Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
% -------------------------
% --- Metodo de Houbolt ---
% -------------------------
% Integrador implicito de 3a ordem (Bathe, Finite Element Procedures).
% Os 3 primeiros instantes (t= 0, dt, 2dt) sao obtidos pela diferenca
% central (MDC2); a partir dai aplica-se a formula de Houbolt.
%
% Variaveis de entrada:
% M - Matriz de massa do sistema
% K - Matriz de rigidez do sistema
% C - Matriz de amortecimento do sistema
% R - Matriz de força aplicada (Ne x Np+1)
% u0 - Deslocamento inicial
% u0p - Velocidade inicial
% u0pp - Aceleracao inicial
% dt - Variacao temporal (Discretizacao)
% Np - Numero de passos
%
% Saidas (Ne x Np+1): u, up, upp nos instantes 0, dt, ..., Np*dt

function [u, up, upp]= MH3(M,K,C,R,u0,u0p,u0pp,dt,Np)
%% Parametros iniciais

Ne= size(u0,1);                             % Numero de elementos

a0= 2/dt^2;                                 % Constante a0
a1= 11/(6*dt);                              % Constante a1
a2= 5/dt^2;                                 % Constante a2
a3= 3/dt;                                   % Constante a3
a4= -2*a0;                                  % Constante a4
a5= -a3/2;                                  % Constante a5
a6= a0/2;                                   % Constante a6
a7= a3/9;                                   % Constante a7

Kc= K+ a0*M+ a1*C;                          % K chapeu

% Partida pela diferenca central: U(:,k+1) -> t= k*dt
[U, UP, UPP]= MDC2(M,K,C,R,u0,u0p,u0pp,Np,dt);

%% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np+1);
up= zeros(Ne, Np+1);
upp= zeros(Ne, Np+1);
Rc= zeros(Ne, Np+1);
u(:,1:3)= U(:,1:3);
up(:,1:3)= UP(:,1:3);
upp(:,1:3)= UPP(:,1:3);

for i= 3:Np
    % R Chapeu
    Rc(:,i+1)= R(:,i+1)+ M*(a2*u(:,i)+a4*u(:,i-1)+a6*u(:,i-2))+ C*(a3*u(:,i)+a5*u(:,i-1)+a7*u(:,i-2));

    % Calculando deslocamentos u_{t+\Delta t}
    u(:,i+1)= Kc\Rc(:,i+1);

    % Calculando acelerações e velocidades - upp_{t+\Delta t} e up_{t+\Delta t}
    up(:,i+1)= a1*u(:,i+1)- a3*u(:,i)- a5*u(:,i-1)- a7*u(:,i-2);
    upp(:,i+1)= a0*u(:,i+1)- a2*u(:,i)- a4*u(:,i-1)- a6*u(:,i-2);

end
end
