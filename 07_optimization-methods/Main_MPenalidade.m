% Main Metodo da Penalidade
clear all; close all; clc;

% Função P1

fx= @(x) 1/2*(x(1)^2 + (x(2)^2)/3);         %Função inicial
g= @(x) [x(1); x(2)/3];                     %Grandiente da função

Pini= @(x) (x(1)+x(2)-1);                   %Penalização
P= @(x) Pini(x)^2;                          %Penalização (quadrática)

dPini= @(x) [1;1];
dP= @(x) 2*Pini(x)*dPini(x);
% dP= @(x) 2*Pini(x) + dPini(x)^2;

x0= [0; 0];                                 %x inicial
xf= [0.25; 0.75];                           %x ótimo final

ci= 2;
gamma= 1.3;
errp= 0.00001;
nim= 1000;
itmaxDico= nim;
results= 1;

x1=0;                                       %Intervalo dicotômica 
x2=1;

[xP,xrvec,i]= MPenalidade(fx,g,P,dP,x0,ci,gamma,x1,x2,nim,errp,itmaxDico,results);

% xfm= fmincon(fx,[0.5, 0.5],[],[],[],[],[0,0],[1,1]);

figure
plot(xrvec(1,:), 'linewidth', 1.5); grid on; grid minor;
title('Convergence of x_1');
xlabel('Iteration'); ylabel('x_1');
ylim([0 0.26])

figure
plot(xrvec(2,:), 'linewidth', 1.5); grid on; grid minor;
title('Convergence of x_2');
xlabel('Iteration'); ylabel('x_2');

% figure
% plot(xrvec(1,:),xrvec(2,:), 'g','linewidth', 1.5); grid on; hold on
% plot(xrvec(1,:),xrvec(2,:), 'bo', 'linewidth', 1.5);
% plot(0.25,0.75, 'rx', 'linewidth', 1.4);
% title('Convergência Final');
% xlabel('Valor de x1'); ylabel('Valor de x2');
% xlim([0 0.26]); ylim([0 0.78]);
% legend('Convergência','Valor da iteração','Mínimo absoluto','Location','SouthEast')

%% Problema 2
clear all; close all; clc;

% Função P1

fx= @(x) exp(x(1)*x(2)*x(3)*x(4)*x(5));         %Função inicial
g1= @(x) x(2)*x(3)*x(4)*x(5)*(exp(x(1)*x(2)*x(3)*x(4)*x(5)));  %Grandiente da função
g2= @(x) x(1)*x(3)*x(4)*x(5)*(exp(x(1)*x(2)*x(3)*x(4)*x(5)));  %Grandiente da função
g3= @(x) x(1)*x(2)*x(4)*x(5)*(exp(x(1)*x(2)*x(3)*x(4)*x(5)));  %Grandiente da função
g4= @(x) x(1)*x(2)*x(3)*x(5)*(exp(x(1)*x(2)*x(3)*x(4)*x(5)));  %Grandiente da função
g5= @(x) x(1)*x(2)*x(3)*x(4)*(exp(x(1)*x(2)*x(3)*x(4)*x(5)));  %Grandiente da função
g= @(x) [g1(x); g2(x); g3(x); g4(x); g5(x)];

Pini1= @(x) (x(1)^2+x(2)^2+x(3)^2+x(4)^2+x(5)^2-10);  %Penalização 1
Pini2= @(x) (x(2)*x(3)-5*x(4)*x(5));        %Penalização 2
Pini3= @(x) (x(1)^3+x(2)^3+1);              %Penalização 3
Pisum= @(x) [Pini1(x), Pini2(x), Pini3(x)];
P= @(x) sum((Pisum(x).*(Pisum(x)~=0)).^2);                          %Penalização (quadrática)

dPini1= @(x) 2*[x(1);x(2);x(3);x(4);x(5)];
dPini2= @(x) [0;x(3);x(2);-5*x(5);-5*x(4)];
dPini3= @(x) 3*[x(1)^2;x(2)^2;0;0;0];
dP= @(x) 2*(Pini1(x)*dPini1(x)+Pini2(x)*dPini2(x)+Pini3(x)*dPini3(x));

x0= [-2 ; 2 ; 2 ; -1 ; -1];
xf= [ -1.717 ; 1.596 ; 1.827 ; -0.764 ; -0.763 ];               %x ótimo final

ci= 2;
gamma= 1.3;
errp= 0.00001;
nim= 1000;
itmaxDico= nim;
results= 0;

x1=0;                                       %Intervalo dicotômica 
x2=1;

[xP,xrvec,i]= MPenalidade(fx,g,P,dP,x0,ci,gamma,x1,x2,nim,errp,itmaxDico,results);

% xfm= fmincon(fx,[0.5, 0.5],[],[],[],[],[0,0],[1,1]);
% 
% figure
% plot(xrvec(1,:), 'linewidth', 1.5); grid on; grid minor;
% title('Convergência x1');
% xlabel('Iteração'); ylabel('Valor de x1');
% ylim([0 0.26])
% 
% figure
% plot(xrvec(2,:), 'linewidth', 1.5); grid on; grid minor;
% title('Convergência x2');
% xlabel('Iteração'); ylabel('Valor de x2');

% figure
% plot(xrvec(1,:),xrvec(2,:), 'g','linewidth', 1.5); grid on; hold on
% plot(xrvec(1,:),xrvec(2,:), 'bo', 'linewidth', 1.5);
% plot(0.25,0.75, 'rx', 'linewidth', 1.4);
% title('Convergência Final');
% xlabel('Valor de x1'); ylabel('Valor de x2');
% xlim([0 0.26]); ylim([0 0.78]);
% legend('Convergência','Valor da iteração','Mínimo absoluto','Location','SouthEast')


