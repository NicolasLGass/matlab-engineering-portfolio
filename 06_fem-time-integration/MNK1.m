%% Elementos Finitos 3 - Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
% ----------------------------
% --- Metodo de Newmark v1 ---
% ----------------------------
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
% delt - Delta escolhido (delt>=0.5)
% alp - Alpha escolhido (alp>=0.25*(0.5+delt)^2)

function [u, up, upp]= MNK1(M,K,C,R,u0,u0p,u0pp,dt,Np,delt,alp)
%% Parametros iniciais

Ne= size(u0,1);                 % Numero de elementos

a0= 1/(alp*dt^2);               % Constante a0
a1= delt/(alp*dt);              % Constante a1
a2= 1/(alp*dt);                 % Constante a2
a3= (1/(2*alp))-1;              % Constante a3
a4= (delt/alp)-1;               % Constante a4
a5= (dt/2)*((delt/alp)-2);      % Constante a5
a6= dt*(1-delt);                % Constante a6
a7= delt*dt;                    % Constante a7

Kc= K+a0*M+a1*C;                % K chapeu

%% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);               % Iniciando os vetores
up= zeros(Ne, Np);
upp= zeros(Ne, Np);
Rc= zeros(Ne, Np);

u(:,1)= u0;                     % Atribuido valores iniciais
up(:,1)= u0p;
upp(:,1)= u0pp;

for i= 1:Np
    % R Chapeu
    Rc(:,i+1)= R(:,i+1)+ M*(a0*u(:,i)+a2*up(:,i)+a3*upp(:,i))+ C*(a1*u(:,i)+a4*up(:,i)+a5*upp(:,i));
    
    % Calculando deslocamentos u_{t+\Delta t} (Mins Quadrados)
    u(:,i+1)= Kc\Rc(:,i+1);
    
    % Calculando acelerações e velocidades upp_{t+\Delta t}, up_{t+\Delta t}
    upp(:,i+1)= a0*(u(:,i+1)-u(:,i))- a2*up(:,i)- a3*upp(:,i);
    up(:,i+1)= up(:,i)+ a6*upp(:,i)+ a7*upp(:,i+1);
    
end
end
%%