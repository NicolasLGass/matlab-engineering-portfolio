%% Elementos Finitos 3 - Trabalho 2 - Diferenca central
% Nicolas Londero Gass
clear all; close all; clc;

%% Carregando dados (Posteriormente transformar em função)

M= [2 0;...                 % Matriz de massa
    0 1];

K= [6 -2;...                % Matriz de rigidez 
    -2 4];

C= [0 0;...                 % Matriz de amortecimento
    0 0];
    
R= [0; 10];                 % Vetor de carregamentos

u0= [0; 0];                 % Deslocamentos iniciais
u0p= [0; 0];                % Velocidades iniciais
u0pp= [0; 10];              % Acelerações iniciais

Tf= 4.45;                   % tempo final 1
dt= 0.28/1;                 % variacao temporal - caso 1
% dt= 28;                     % variacao temporal - caso 2
Np= 11;

%% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 1/dt^2;                 % Constante a1
a1= (1/2)*1/dt;             % Constante a2
a2= 2*a0;                   % Constante a3
a3= 1/a2;                   % Constante a4

udt= u0-dt*u0p+a3*u0pp;     % u(-dt)
Mc= a0*M+ a1*C;             % M chapeu

% Mc= inv(Mc);               % Matriz M chapeu triangularizada (L*D*L^t)
%certo???

%% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);           % Inicializando o vetor de deslocamento
Rc= zeros(Ne, Np);          % Inicializando o vetor de carregamento
u(:,1)= udt(:);
u(:,2)= u0;

for i= 2:Np+2

    Rc(:,i)= R(:,1)-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    u(:,i+1)= Mc\Rc(:,i);   % Calculando o deslocamento em t+(delta)t

end


