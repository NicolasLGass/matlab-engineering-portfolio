%% Elementos Finitos 3 - Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
% -----------------------------------
% --- Metodo de Wilson theta θ v1 ---
% -----------------------------------
% Variaveis de entrada:
% M - Matriz de massa do sistema
% K - Matriz de rigidez do sistema
% C - Matriz de amortecimento do sistema
% R - Matriz de força aplicada
% u0 - Deslocamento inicial
% u0p - Velocidade inicial
% u0pp - Aceleracao inicial
% dt - Variacao temporal (Discretizacao)
% Np - Numero de pontos
% tht - Valor de theta

function [u, up, upp]= MWT1(M,K,C,R,u0,u0p,u0pp,dt,Np,tht)
%% Parametros iniciais

% Np=Np+1;
Ne= size(u0,1);             % Numero de elementos

a0= 6/(tht*dt)^2;           % Constante a0
a1= 3/(tht*dt);             % Constante a1
a2= 2*a1;                   % Constante a2
a3= tht*dt/2;               % Constante a3
a4= a0/tht;                 % Constante a4
a5= -a2/tht;                % Constante a5
a6= 1-3/tht;                % Constante a6
a7= dt/2;                   % Constante a7
a8= (dt^2)/6;               % Constante a8

Kc= K+a0*M+a1*C;            % K chapeu

%% Calculo dos carregamentos Rc e deslocamento u
Np= Np+1;
Rctht= zeros(Ne, Np);       % Inicializando vetores
u= zeros(Ne, Np);
up= zeros(Ne, Np);
upp= zeros(Ne, Np);
utht= zeros(Ne, Np+1);

u(:,1)= u0;
up(:,1)= u0p;
upp(:,1)= u0pp;

for i= 1:Np-1
    % R Chapeu Theta
    Rctht(:,i)= R(:,i)+ tht*(R(:,i+1)-R(:,i))+ M*(a0*u(:,i)+a2*up(:,i)+2*upp(:,i))+ C*(a1*u(:,i)+2*up(:,i)+a3*upp(:,i));
    
    % Calculando deslocamentos u_{t+\Delta t} (Mins Quadrados)
    utht(:,i+1)= Kc\Rctht(:,i);
        
    % Calculando acelerações e velocidades e deslocamentos - upp_{t+\Delta t}, up_{t+\Delta t} e u_{t+\Delta t}
    upp(:,i+1)= a4*(utht(:,i+1)-u(:,i))+ a5*up(:,i)+ a6*upp(:,i);    
    up(:,i+1)= up(:,i)+ a7*(upp(:,i+1)+upp(:,i));    
    u(:,i+1)= u(:,i)+ dt*up(:,i)+ a8*(upp(:,i+1)+ 2*upp(:,i));

end
end
%%