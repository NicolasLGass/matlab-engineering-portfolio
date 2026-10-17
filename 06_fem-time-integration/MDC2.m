%% Elementos Finitos 3 - Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
% --------------------------------------
% --- Metodo da diferenca central v2 ---
% --------------------------------------
% Variaveis de entrada:
% M - Matriz de massa do sistema
% K - Matriz de rigidez do sistema
% C - Matriz de amortecimento do sistema
% R - Matriz de força aplicada
% u0 - Deslocamento inicial
% u0p - Velocidade inicial
% u0pp - Aceleracao inicial
% Tf - Tempo final
% dt - Variacao temporal (Discretizacao)

function [uf, up, upp]= MDC2(M,K,C,R,u0,u0p,u0pp,Np,dt)

Np= Np+2;                   % Numero de passos

%% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 1/dt^2;                 % Constante a1
a1= (1/2)*1/dt;             % Constante a2
a2= 2*a0;                   % Constante a3
a3= 1/a2;                   % Constante a4

udt= u0-dt*u0p+a3*u0pp;     % u(-dt)
Mc= a0*M + a1*C;            % M chapeu

%% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);
up= zeros(Ne, Np);
upp= zeros(Ne, Np);
Rc= zeros(Ne, Np);
u(:,1)= udt(:);
u(:,2)= u0;

% Coluna i de u corresponde a t= (i-2)*dt -> carga R(:,i-1)
for i= 2:Np
    % R Chapeu
    Rc(:,i)= R(:,min(i-1,end))-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    
    % Calculando deslocamentos u_{t+\Delta t} (Mins Quadrados)
    u(:,i+1)= Mc\Rc(:,i);
    
    % Calculando acelerações e velocidades - upp_{t+\Delta t} e up_{t+\Delta t}
    upp(:,i)= a0*(u(:,i-1)-2*(u(:,i))+u(:,i+1));
    up(:,i)= a1*(-u(:,i-1)+u(:,i+1));

end

% Saidas alinhadas nos instantes 0, dt, ..., Np*dt
uf= u(:,2:Np);
up= up(:,2:Np);
upp= upp(:,2:Np);
end
