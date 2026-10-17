%% Metodo da Diferenca central
% Variaveis de entrada:
% M - Matriz de massa do sistema
% K - Matriz de rigidez do sistema
% C - Matriz de amortecimento do sistema
% R - Matriz de força aplicada
% u0 - Deslocamento inicial
% u0p - Velocidade inicial
% u0pp - Aceleracao inicial
% Tf - Tempo final

function u= MDC1(M,K,C,R,u0,u0p,u0pp,Tf,dt)

Np= round(Tf/dt+1);         % Numero de passos (int8 saturava em 127)

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

u= zeros(Ne, Np);
Rc= zeros(Ne, Np);
u(:,1)= udt(:);
u(:,2)= u0;

for i= 2:Np+2
    
    Rc(:,i)= R(:,1)-(K-a2*M)*u(:,i) - (a0*M-a1*C)*u(:,i-1);
    u(:,i+1)= Mc\Rc(:,i);

end

% u

