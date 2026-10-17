%% Elementos Finitos 3 - Análise de estabilidade dos métodos de integração
% Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
clear all; close all; clc;
warning off;

%% Carregando dados

Tf= 2.8;                        % Periodo critico
dt= .28;                        % Variacao temporal - caso 1
Np= 12;                         % Numero de pontos
dtcd= dt/pi;                    % Delta t critico (Dif Central)

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

% Deslocamento, Velocidade, e Aceleracao inicias
u0= [0; 0];
u0p= [0; 0];
u0pp= [0; 10];

% Variaveis dos integradores
tht= 1.4;                       % Theta
delt= 0.5;                      % Delta (delt>=0.5)
alp= 0.25;                      % Alpha (alp>=0.25*(0.5+delt)^2)

[un1, upn1, uppn1]= MNK1(M,K,C,R,u0,u0p,u0pp,dt,Np,delt,alp);       % Metodo de Newmark
[uw1, upw1, uppw1]= MWT1(M,K,C,R,u0,u0p,u0pp,dt,Np,tht);            % Metodo de Wilson θ
[uh1, uph1, upph1]= MH3(M,K,C,R,u0,u0p,u0pp,dt,Np);              % Metodo de Houbolt
[ud1, upd1, uppd1]= MDC2(M,K,C,R,u0,u0p,u0pp,Np,dt);                % Metodo de Diferenca Central
                        
%% 1 Modificando dt 

dt=Tf/(pi*6);

[un2, upn2, uppn2]= MNK1(M,K,C,R,u0,u0p,u0pp,dt,Np,delt,alp);       % Metodo de Newmark
[uw2, upw2, uppw2]= MWT1(M,K,C,R,u0,u0p,u0pp,dt,Np,tht);            % Metodo de Wilson θ
[uh2, uph2, upph2]= MH3(M,K,C,R,u0,u0p,u0pp,dt,Np);              % Metodo de Houbolt
[ud2, upd2, uppd2]= MDC2(M,K,C,R,u0,u0p,u0pp,Np,dt);                % Metodo de Diferenca Central
   
%% 2 Modificando dt 

dt=Tf/(pi/1);

[un3, upn3, uppn3]= MNK1(M,K,C,R,u0,u0p,u0pp,dt,Np,delt,alp);       % Metodo de Newmark
[uw3, upw3, uppw3]= MWT1(M,K,C,R,u0,u0p,u0pp,dt,Np,tht);            % Metodo de Wilson θ
[uh3, uph3, upph3]= MH3(M,K,C,R,u0,u0p,u0pp,dt,Np);              % Metodo de Houbolt
[ud3, upd3, uppd3]= MDC2(M,K,C,R,u0,u0p,u0pp,Np,dt);                % Metodo de Diferenca Central
   
%% 3 Modificando dt 

dt=Tf/(pi/6);

[un4, upn4, uppn4]= MNK1(M,K,C,R,u0,u0p,u0pp,dt,Np,delt,alp);       % Metodo de Newmark
[uw4, upw4, uppw4]= MWT1(M,K,C,R,u0,u0p,u0pp,dt,Np,tht);            % Metodo de Wilson θ
[uh4, uph4, upph4]= MH3(M,K,C,R,u0,u0p,u0pp,dt,Np);              % Metodo de Houbolt
[ud4, upd4, uppd4]= MDC2(M,K,C,R,u0,u0p,u0pp,Np,dt);                % Metodo de Diferenca Central

%% Comp Metodos u_{t1} - dt1
ax= 0:Np;

figure(1)
plot(ax, ud1(1,:)); hold on; grid on;
plot(ax, uh1(1,:));
plot(ax, uw1(1,:));
plot(ax, un1(1,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'northwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t1} (\Deltat = 0.28)', 'fontsize', 15)

% u_{t2} - dt1
figure(2)
plot(ax, ud1(2,:)); hold on; grid on;
plot(ax, uh1(2,:));
plot(ax, uw1(2,:));
plot(ax, un1(2,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'northwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t2} (\Deltat = 0.28)', 'fontsize', 15)

%% Comp Metodos u_{t1} - dt2
ax= 0:Np;

figure(3)
plot(ax, ud2(1,:)); hold on; grid on;
plot(ax, uh2(1,:));
plot(ax, uw2(1,:));
plot(ax, un2(1,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'northwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t1} (\Deltat = T_f/(\pi*6))', 'fontsize', 15)

% u_{t2} - dt2
figure(4)
plot(ax, ud2(2,:)); hold on; grid on;
plot(ax, uh2(2,:));
plot(ax, uw2(2,:));
plot(ax, un2(2,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'northwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t2} (\Deltat = T_f/(\pi*3))', 'fontsize', 15)

%% Comp Metodos u_{t1} - dt3
ax= 0:Np;

figure(5)
plot(ax, ud3(1,:)); hold on; grid on;
plot(ax, uh3(1,:));
plot(ax, uw3(1,:));
plot(ax, un3(1,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'northwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t1} (\Deltat = T_f/(\pi))', 'fontsize', 15)

% u_{t2} - dt3
figure(6)
plot(ax, ud3(2,:)); hold on; grid on;
plot(ax, uh3(2,:));
plot(ax, uw3(2,:));
plot(ax, un3(2,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'northwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t2} (\Deltat = T_f/(\pi))', 'fontsize', 15)

%% Comp Metodos u_{t1} - dt4
ax= 0:Np;

figure(7)
plot(ax, ud4(1,:)); hold on; grid on;
plot(ax, uh4(1,:)); hold on; grid on;
plot(ax, uw4(1,:));
plot(ax, un4(1,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'northwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t1} (\Deltat = T_f/(\pi/6))', 'fontsize', 15)

% u_{t2} - dt4
figure(8)
plot(ax, ud4(2,:)); hold on; grid on;
plot(ax, uh4(2,:)); hold on; grid on;
plot(ax, uw4(2,:));
plot(ax, un4(2,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'southwest', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t2} (\Deltat = T_f/(\pi/6))', 'fontsize', 15)

%% Comp Metodos u_{t1} - dt4 zoom
ax= 0:Np;

figure(9)
plot(ax, ud4(1,:)); hold on; grid on;
plot(ax, uh4(1,:)); hold on; grid on;
plot(ax, uw4(1,:), 'g');
plot(ax, un4(1,:));
legend('Central diff.','Houbolt','Wilson \theta','Newmark', 'location', 'southeast', 'fontsize', 10);
xlabel('Time step n', 'fontsize', 13); ylabel('Displacement amplitude', 'fontsize', 13);
title('Comparison of results u_{t1} (\Deltat = T_f/(\pi/6))', 'fontsize', 15)
axis([0 12 -3 5])
    
% u_{t2} - dt4
% figure(10)
% plot(ax, ud4(2,:)); hold on; grid on;
% plot(ax, uh4(2,:)); hold on; grid on;
% plot(ax, uw4(2,:));
% plot(ax, un4(2,:));
% legend('Dif Central','Houbolt','Wilson \theta','Newmark', 'location', 'southwest', 'fontsize', 10);
% xlabel('Passo \deltat', 'fontsize', 13); ylabel('Amplitude de deslocamento', 'fontsize', 13);
% title('Comparação deresultados u_{t2} (\Deltat = T_f/(\pi/6))', 'fontsize', 15)
% axis([-5 5 0 12])

%%
% %% Newmark
% ax= 0:Np;
% 
% figure(1)
% plot(ax, un1(1,:)); hold on; grid on;
% plot(ax, un1(2,:));
% 
% figure(2)
% plot(ax, un2(1,:)); hold on; grid on;
% plot(ax, un2(2,:));
% 
% figure(3)
% plot(ax, un3(1,:)); hold on; grid on;
% plot(ax, un3(2,:));
% 
% figure(4)
% plot(ax, un4(1,:)); hold on; grid on;
% plot(ax, un4(2,:));
% 
% %% Houbolt
% ax= 0:Np;
% 
% figure(5)
% plot(ax, uh1(1,:)); hold on; grid on;
% plot(ax, uh1(2,:));
% 
% figure(6)
% plot(ax, uh2(1,:)); hold on; grid on;
% plot(ax, uh2(2,:));
% 
% figure(7)
% plot(ax, uh3(1,:)); hold on; grid on;
% plot(ax, uh3(2,:));
% 
% figure(8)
% plot(ax, uh4(1,:)); hold on; grid on;
% plot(ax, uh4(2,:));
% 
% %% Wilson
% ax= 0:Np;
% 
% figure(9)
% plot(ax, uw1(1,:)); hold on; grid on;
% plot(ax, uw1(2,:));
% 
% figure(10)
% plot(ax, uw2(1,:)); hold on; grid on;
% plot(ax, uw2(2,:));
% 
% figure(11)
% plot(ax, uw3(1,:)); hold on; grid on;
% plot(ax, uw3(2,:));
% 
% figure(12)
% plot(ax, uw4(1,:)); hold on; grid on;
% plot(ax, uw4(2,:));
% 
% %% Dif Central
% ax= 0:Np;
% 
% figure(13)
% plot(ax, ud1(1,:)); hold on; grid on;
% plot(ax, ud1(2,:));
% 
% figure(14)
% plot(ax, ud2(1,:)); hold on; grid on;
% plot(ax, ud2(2,:));
% 
% figure(15)
% plot(ax, ud3(1,:)); hold on; grid on;
% plot(ax, ud3(2,:));
% 
% figure(16)
% plot(ax, ud4(1,:)); hold on; grid on;
% plot(ax, ud4(2,:));
