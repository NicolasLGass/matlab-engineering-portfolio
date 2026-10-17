%% Elementos Finitos 3 - Trabalho 2 - Diferenca central
% Nicolas Londero Gass
clear all; close all; clc;

% Carregando dados (Posteriormente transformar em função)

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

T1= 4.45;                   % Período freq 1
T2= 2.8;                    % Período freq 2 (crítica, por ser menor/freq mais alta)
dt= 0.28/1;                 % variacao temporal - caso 1
Np= 11;

% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 1/dt^2;                 % Constante a1
a1= (1/2)*1/dt;             % Constante a2
a2= 2*a0;                   % Constante a3
a3= 1/a2;                   % Constante a4

udt= u0-dt*u0p+a3*u0pp;     % u(-dt)
Mc= a0*M+ a1*C;             % M chapeu

% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);           % Inicializando o vetor de deslocamento
Rc= zeros(Ne, Np);          % Inicializando o vetor de carregamento
u(:,1)= udt(:);
u(:,2)= u0;

for i= 2:Np+2
    
    Rc(:,i)= R(:,1)-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    u(:,i+1)= Mc\Rc(:,i);   % Calculando o deslocamento em t+(delta)t

end

disp('Valor calculado com dt= 2.8/10 (Estavel): \n')
u

%% Delta t critico (Estavel)

dtc= T2/pi;                 % delta t crítico
dt1= dtc;                   % variacao temporal - estabilidade
Np= 11;

% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 1/dt1^2;                 % Constante a1
a1= (1/2)*1/dt1;             % Constante a2
a2= 2*a0;                   % Constante a3
a3= 1/a2;                   % Constante a4

udt= u0-dt1*u0p+a3*u0pp;     % u(-dt)
Mc= a0*M+ a1*C;             % M chapeu

% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);           % Inicializando o vetor de deslocamento
Rc= zeros(Ne, Np);          % Inicializando o vetor de carregamento
u(:,1)= udt(:);
u(:,2)= u0;

for i= 2:Np+2

    Rc(:,i)= R(:,1)-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    u(:,i+1)= Mc\Rc(:,i);   % Calculando o deslocamento em t+(delta)t

end

disp('Valor calculado com dt= dtc (critico)(Estavel): \n')
uc=u

%% Delta t critico - (2.8/10) (Estavel)

dtc= T2/pi;                 % delta t crítico
dtc2= dtc-T2/10;            % delta t crítico - 0.28
dt= dtc2;                   % variacao temporal - estabilidade
Np= 11;

% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 1/dt^2;                 % Constante a1
a1= (1/2)*1/dt;             % Constante a2
a2= 2*a0;                   % Constante a3
a3= 1/a2;                   % Constante a4

udt= u0-dt*u0p+a3*u0pp;     % u(-dt)
Mc= a0*M+ a1*C;             % M chapeu

% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);           % Inicializando o vetor de deslocamento
Rc= zeros(Ne, Np);          % Inicializando o vetor de carregamento
u(:,1)= udt(:);
u(:,2)= u0;

for i= 2:Np+2

    Rc(:,i)= R(:,1)-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    u(:,i+1)= Mc\Rc(:,i);   % Calculando o deslocamento em t+(delta)t

end

disp('Valor calculado com dt= dtc-(2.8/10) (Estavel): \n')
uc1=u

%% Delta t crítico + (2.8/10) (Instavel)

dtc3= dtc+T2/10;            % delta t crítico + 0.28
dt= dtc3;                   % variacao temporal - estabilidade
Np= 11;

% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 1/dt^2;                 % Constante a1
a1= (1/2)*1/dt;             % Constante a2
a2= 2*a0;                   % Constante a3
a3= 1/a2;                   % Constante a4

udt= u0-dt*u0p+a3*u0pp;     % u(-dt)
Mc= a0*M+ a1*C;             % M chapeu

% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);           % Inicializando o vetor de deslocamento
Rc= zeros(Ne, Np);          % Inicializando o vetor de carregamento
u(:,1)= udt(:);
u(:,2)= u0;

for i= 2:Np+2

    Rc(:,i)= R(:,1)-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    u(:,i+1)= Mc\Rc(:,i);   % Calculando o deslocamento em t+(delta)t

end

disp('Valor calculado com dt= dtc+(2.8/10) (Instavel): \n')
uc2=u

%% Delta t crítico + vmi (vmc = Valor Minimo de acrescimo para Instabilidade - Estimado por tentativa e erro)

dtvmi= dtc+0.0450398;                   % variacao temporal - estabilidade
Np= 11;

% Parametros iniciais

Ne= size(u0,1);             % Numero de elementos

a0= 1/dtvmi^2;                 % Constante a1
a1= (1/2)*1/dtvmi;             % Constante a2
a2= 2*a0;                   % Constante a3
a3= 1/a2;                   % Constante a4

udt= u0-dtvmi*u0p+a3*u0pp;     % u(-dt)
Mc= a0*M+ a1*C;             % M chapeu

% Calculo dos carregamentos Rc e deslocamento u

u= zeros(Ne, Np);           % Inicializando o vetor de deslocamento
Rc= zeros(Ne, Np);          % Inicializando o vetor de carregamento
u(:,1)= udt(:);
u(:,2)= u0;

for i= 2:Np+2

    Rc(:,i)= R(:,1)-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    u(:,i+1)= Mc\Rc(:,i);   % Calculando o deslocamento em t+(delta)t

end

disp('Valor calculado com dt= dtc+vmi (vmc = Valor Minimo de acrescimo para Instabilidade - Estimado por tentativa e erro): \n')
uc3=u


