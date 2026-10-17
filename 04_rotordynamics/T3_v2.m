% Analise de maquinas rotativas - Exercicio 3 - Aproximacao Modal 3GDL
% Nicolas Londero Gass
clear all; close all; clc
if exist('OCTAVE_VERSION', 'builtin'), pkg load symbolic; end   % Octave: symbolic package
warning off

% Sistema do problema
%        >> x1(t)    >> x2(t)    >> x3(t)   
%      ___|__      ___|__      ___|__     
% -K1-|      |-K2-|      |-K3-|      |-K4-
%     |  M1  |    |  M2  |    |  M3  |    
% -C1-|______|-C2-|______|-C3-|______|-C4-

%% Dados do sistema de 3GDL

M= [2 0 0;...                       % Matriz de massa [kg]
    0 3 0;...
    0 0 2];

C= [4.2  -2    0;...                % Matriz de amortecimento [N.s/m]
    -2   6.3  -4;...
     0   -4   8.2];

K= (10^5)*[4 -2 0;...               % Matriz de rigidez [N/m]
          -2 6 -4;...
           0 -4 8];

T= [0 0 0;...                       % Matriz temporária zero
    0 0 0;...
    0 0 0];

%% Calculando os autovetores e autovalores

CM= [C M;...                        % Matriz C e M
     M T];
 
KM= [K T;...                        % Matriz K e M                      
     T -M];

% form => formas modais (phi) freqm => freqs modais ()
[form, freqm]= eig(KM, -CM);        % calculando frequencias modais(freqm) e formas modais (form)

form= form(4:6, 5:-2:1);            % Selecionando frequencias modais
phi= real(form);                    % Matriz de formas modais

freqm2= diag(freqm);                % Selecionando autovalor passo 1
freqm2= freqm2(5:-2:1);             % Selecionando autovalor passo 2 (Si= sigma_i + j*omgd_i)

omgd= diag(imag(freqm2));           % Frequencia natural amortecida (omgd_i= omgn_i*sqrt(1-qsi_i^2))
omgdv= imag(freqm2);                % Vetor de frequencias naturais amortecidas
sigma= real(freqm2);                % Vetor de sigma (sigma= qsi_i*omgn_i)

for i=1:length(freqm2)
    
    fatC(i)= sqrt(1/(1+((omgdv(i)/sigma(i))^2)));       % Fator de amortecimento (qsi)
    omgn(i)= abs(sigma(i)/fatC(i));                     % Frequencias naturais
                                                                    
    gamm(i)= 1/sqrt(phi(:,i)'*M*phi(:,i));              % Variavel gamma (Vetor modal)
    Phin(:,i)= gamm(i)*phi(:,i);                        % Phi normalizado
end

sigmaComp= -(fatC.*omgn)';          % Verificacao do valor de sigma por outro calculo

if sigmaComp-sigma<=1e-8            % Verificando convergencia de sigma
    aSigConv= 'Valores iguais';
else
    aSigConv= 'Valores diferentes';
end

% Normalizando as matrizes do sistema

Wn2= Phin'*K*Phin;                  % Matriz de freqs naturais (wn^2)
Cn= Phin'*C*Phin;                   % Matriz de amortecimento normalizada
Mn= Phin'*M*Phin;                   % Matriz de massa normalizada (I)

%% Calculando o fator de amortecimento 'h' (qsi, deve dar igual a variavel fatC)

omg= sqrt(diag(Wn2));               % Freqs naturais (verificando valores)
ceqc= 2*omg;                        % Fator de amortecimento equivalente critico (valido apenas para massa modal normalizada = 1) 
ceq= diag(Cn);                      % Fator de amortecimento equivalente
h= (ceq./ceqc);                     % Fator de amortecimento (qsi)

if h-fatC'<=1e-8                    % Verificando convergencia de qsi
    aQsiConv= 'Valores iguais';
else
    aQsiConv= 'Valores diferentes';
end

%% Alpha e Beta por minimos quadrados - [C]= alp*[M]+bet*[K]

syms alp bet
E= (h-(alp/2)./omg)-((bet/2).*omg);	% Função Erro (usar 'h' ou 'fatC')
Es= sum(E.^2);

eqalp= diff(Es, alp)==0;
eqbet= diff(Es, bet)==0;
[alpt, bett]= solve([eqalp, eqbet], [alp bet]);
alpf= double(alpt);                   % Alpha final
betf= double(bett);                   % Beta final

Cf= alpf*M+betf*K;                  % Matriz de amortecimento C final 
Cf2= alpf*Mn+betf*Wn2;              % Matriz de amortecimento Cn

if (C-Cf)<=1e-8                     % Verificando convergencia da matriz C
    aCConv= 'Valores iguais';
else
    aCConv= 'Valores diferentes';
end

if (Cn-Cf2)<=1e-8                   % Verificando convergencia da matriz Cn
    aCnConv= 'Valores iguais';
else
    aCnConv= 'Valores diferentes';
end

%% Aproximacao modal - x(t)=Xn.q(t)

tvec= 0:0.0001:10;                  % Vetor de tempo
q= (Phin^-1)*[0; 0; 1];             % Deslocamento modal
qp= (Phin^-1)*[0; 0; 0];            % Aceleracao modal
cc= [q; qp];                        % Cond Iniciais (vetor coluna): x1, x2, x3, xp1, xp2, xp3

[t, sol]= ode45(@(tvec,q)eqsp(Wn2, Cn, q, tvec),tvec,cc);

qf= sol(:,1:3);
qpf= sol(:,4:6);

xt= Phin*qf';                       % x(t) em coordenadas fisicas

%% Resultados coordenadas fisica
% fname = fullfile(pwd, 'output_figures');   % Descomente (junto com os saveas) para exportar as figuras
fontsize= 29;

figure(1)
plot(t,xt(1,:), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Displacement x_1 - physical coordinates')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'x1'), 'png');

figure(2)
plot(t,xt(2,:), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Displacement x_2 - physical coordinates')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'x2'), 'png');

figure(3)
plot(t,xt(3,:), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Displacement x_3 - physical coordinates')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'x3'), 'png');

figure(4)
plot(t,(xt(1,:)+xt(2,:)+xt(3,:)), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Sum of physical displacements x')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'xs'), 'png');

%% Resultados coordenadas modais

figure(5)
plot(t,qf(:,1), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Displacement q_1 - modal coordinates')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'q1'), 'png');

figure(6)
plot(t,qf(:,2), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Displacement q_2 - modal coordinates')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'q2'), 'png');

figure(7)
plot(t,qf(:,3), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Displacement q_3 - modal coordinates')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'q3'), 'png');

figure(8)
plot(t,(qf(:,1)+qf(:,2)+qf(:,3)), 'linewidth', 0.8); grid on; grid minor
xlabel('Time [s]'); ylabel('Amplitude')
title('Sum of modal displacements q')
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'FontSize',fontsize)
% saveas(gca, fullfile(fname, 'qs'), 'png');

%%