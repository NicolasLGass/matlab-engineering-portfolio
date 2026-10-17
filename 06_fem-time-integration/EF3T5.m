%% Elementos Finitos 3 - Trabalho 5 - Metodo Newmark
% Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
clear all; close all; clc;

%% Carregando dados (Posteriormente transformar em funcao)

Tf= 2.8;                        % Tempo final
dt= .28;                        % Variacao temporal - caso 1
Np= 12;                         % Numero de pontos

M= [2 0;...                     % Matriz de massa
    0 1];

K= [6 -2;...                    % Matriz de rigidez
    -2 4];

C= [0 0;...                     % Matriz de amortecimento
    0 0];

rt= [0; 10];
for i=1:Np+1
    R(:,i)= rt;                 % Matriz de forcas (p/ cada intervalo de tempo)
end

u0= [0; 0];
u0p= [0; 0];
u0pp= [0; 10];

delt= 0.5;                      % Delta (delt>=0.5)
alp= 0.25;                      % Alpha (alp>=0.25*(0.5+delt)^2)
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
    % R chapeu
    Rc(:,i+1)= R(:,i+1)+ M*(a0*u(:,i)+a2*up(:,i)+a3*upp(:,i))+ C*(a1*u(:,i)+a4*up(:,i)+a5*upp(:,i));
    
    % Calculando deslocamentos u_{t+\Delta t} (Mins Quadrados)
    u(:,i+1)= Kc\Rc(:,i+1);
    
    % Calculando acelerações e velocidades upp_{t+\Delta t}, up_{t+\Delta t}
    upp(:,i+1)= a0*(u(:,i+1)-u(:,i))- a2*up(:,i)- a3*upp(:,i);
    up(:,i+1)= up(:,i)+ a6*upp(:,i)+ a7*upp(:,i+1);
    
end

u

%%