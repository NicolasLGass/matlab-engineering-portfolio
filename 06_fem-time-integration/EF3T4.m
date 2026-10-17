%% Elementos Finitos 3 - Trabalho 4 - Metodo de Wilson theta θ
% Nicolas Londero Gass
clear all; close all; clc;

%% Carregando dados (Posteriormente transformar em funcao)

Tf= 2.8;            % Tempo final
dt= .28;            % Variacao temporal - caso 1
Np= 13;             % Numero de pontos
tht= 1.4;           % Theta

M= [2 0;...         % Matriz de massa
    0 1];

K= [6 -2;...        % Matriz de rigidez
    -2 4];

C= [0 0;...         % Matriz de amortecimento
    0 0];

rt= [0; 10];
for i=1:Np+1
    R(:,i)= rt;    % Matriz de forcas (p/ cada intervalo de tempo)
end

u0= [0; 0];
u0p= [0; 0];
u0pp= [0; 10];

%% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 6/(tht*dt)^2;           % Constante a0
a1= 3/(tht*dt);             % Constante a1
a2= 2*a1;                   % Constante a2
a3= tht*dt/2;             % Constante a3
a4= a0/tht;                 % Constante a4
a5= -a2/tht;                % Constante a5
a6= 1-3/tht;                % Constante a6
a7= dt/2;                   % Constante a7
a8= (dt^2)/6;               % Constante a8

Kc= K+a0*M+a1*C;            % K chapeu

%% Calculo dos carregamentos Rc e deslocamento u
Np= Np+0;
Rctht= zeros(Ne, Np);          % Inicializando vetores
u= zeros(Ne, Np);
up= zeros(Ne, Np);
upp= zeros(Ne, Np);

u(:,1)= u0;
up(:,1)= u0p;
upp(:,1)= u0pp;

for i= 1:Np-1
    
    Rctht(:,i)= R(:,i)+ tht*(R(:,i+1)-R(:,i))+ M*(a0*u(:,i)+a2*up(:,i)+2*upp(:,i))+ C*(a1*u(:,i)+2*up(:,i)+a3*upp(:,i));          % R chapeu
    
    utht(:,i+1)= Kc\Rctht(:,i);
        
    upp(:,i+1)= a4*(utht(:,i+1)-u(:,i))+ a5*up(:,i)+ a6*upp(:,i);
    
    up(:,i+1)= up(:,i)+ a7*(upp(:,i+1)+upp(:,i));
    
    u(:,i+1)= u(:,i)+ dt*up(:,i)+ a8*(upp(:,i+1)+ 2*upp(:,i));

end

u
up
upp

%%