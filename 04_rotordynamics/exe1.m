% Analise de maquinas rotativas - Exercicio 1 
% Nicolas Londero Gass
clear all; close all; clc

%% Solução da equacoes do movimento - Primeiro caso
t= 0:0.0001:2;                      %Vetor temporal
CC= [0;0;0;0;0;0];                  %Condições de contorno
[t, res]= ode45(@EDT1, t, CC);      %Solver da equação
res(:,6)= (60/(2*pi))*res(:,6);    %rad/s para rpm

%% Gráfico 1 - Torque = 1.834
figure(1)
subplot(2,1,1);
plot_disp_speed(t, res(:,1), res(:,6), [-0.0055 0.0055], 'Displacement u_y(t) [m]', 'Displacements u_y and u_z - torque 1.834 N·m');

subplot(2,1,2);
plot_disp_speed(t, res(:,2), res(:,6), [-0.0055 0.0055], 'Displacement u_z(t) [m]');

MaxC11= max(abs(res(:,1)));
MaxC12= max(abs(res(:,2)));
MaxT1= max(res(:,6));
%%
%% Solução da equacoes do movimento - Segundo caso
t= 0:0.0001:5;                      %Vetor temporal
CC= [0;0;0;0;0;0];                  %Condições de contorno
[t, res]= ode45(@EDT2, t, CC);      %Solver da equação
res(:,6)= (60/(2*pi))*res(:,6);    %%rad/s para rpm

%% Gráfico 2 - Torque = 0.0825
figure(2)
subplot(2,1,1);
plot_disp_speed(t, res(:,1), res(:,6), [-0.012 0.012], 'Displacement u_y(t) [m]', 'Displacements u_y and u_z - torque 0.0825 N·m');

subplot(2,1,2);
plot_disp_speed(t, res(:,2), res(:,6), [-0.012 0.012], 'Displacement u_z(t) [m]');

MaxC21= max(abs(res(:,1)));
MaxC22= max(abs(res(:,2)));
MaxT2= max(res(:,6));
%%
%% Solução da equacoes do movimento - Terceiro caso
t= 0:0.0001:5;                      %Vetor temporal
CC= [0;0;0;0;0;0];                  %Condições de contorno
[t, res]= ode45(@EDT3, t, CC);      %Solver da equação
res(:,6)= (60/(2*pi))*res(:,6);    %rad/s para rpm

%% Gráfico 3 - Torque = 0.075
figure(3)
subplot(2,1,1);
plot_disp_speed(t, res(:,1), res(:,6), [-0.014 0.014], 'Displacement u_y(t) [m]', 'Displacements u_y and u_z - torque 0.075 N·m');

subplot(2,1,2);
plot_disp_speed(t, res(:,2), res(:,6), [-0.014 0.014], 'Displacement u_z(t) [m]');

MaxC31= max(abs(res(:,1)));
MaxC32= max(abs(res(:,2)));
MaxT3= max(res(:,6));

%% Comparando resultados entre os casos
fprintf('Máximo deslocamento em y: \n');
fprintf('Caso 1: %.06f m | Caso 2: %.06f m | Caso 3: %.06f m\n\n', MaxC11, MaxC21, MaxC31);
fprintf('Máximo deslocamento em z: \n');
fprintf('Caso 1: %.06f m | Caso 2: %.06f m | Caso 3: %.06f m\n\n', MaxC12, MaxC22, MaxC32);
fprintf('Rotação máxima atingida: \n');
fprintf('Caso 1: %.02f rpm | Caso 2: %.02f rpm | Caso 3: %.02f rpm \n', MaxT1, MaxT2, MaxT3);

