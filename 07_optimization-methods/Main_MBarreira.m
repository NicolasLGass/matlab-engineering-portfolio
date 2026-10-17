% Main Metodo da Penalidade
clear all; close all; clc; format long;

% Função P1

fx= @(x) (x(1)-2)^4 + (x(1)-2*x(2))^2;      %Função inicial
g= @(x) [(4*(x(1)-2)^3) + 2*(x(1)-2*x(2)); -4*(x(1)-2*x(2))];   %Grandiente da função

Bini= @(x) (x(1)^2-x(2));
B= @(x) (1/(Bini(x)));                        %Barreira
% B= @(x) log(-g(x));

dB= @(x) [-(2*x(1))/(x(1)^2-x(2))^2; 1/(x(1)^2-x(2))^2];

x0= [1; 4];                                 %x inicial
% xf= [0.25; 0.75];                         %x ótimo final

ci= 10;                                     %Ci inicial
gamma= 10;
errp= 0.00000000001;                           %Erro de tolerância
nim= 100;                                  %Número máx de iterações
itmaxDico= nim;                             %N max de iterações dicotômica
results= 0;                                 %Mostrar resultados, 1=sim

x1=0;                                       %Intervalo dicotômica 
x2=10;

[xP,xrvec,i,vecval, vecalfa]= MBarreira(fx,g,B,dB,x0,ci,gamma,x1,x2,nim,errp,itmaxDico,results);

figure
plot(xrvec(1,:), 'k-x', 'linewidth', 1.5); grid on; grid minor; hold on;
title('Convergence of x_1');
plot(xrvec(2,:), 'b', 'linewidth', 1.5); grid on; grid minor; hold off;
title('Convergence of x_1 and x_2');
legend('x1','x2')
xlabel('Iteration'); ylabel('x_1, x_2');

xaux= -8:0.00001:8;
for i=1:length(xaux)
    fxp(i)= (xaux(i).^2);
end

figure
plot(xaux,fxp, 'r', 'linewidth', 1); grid on; hold on
plot(xrvec(1,:),xrvec(2,:), 'g','linewidth', 1.5); 
plot(xrvec(1,:),xrvec(2,:), 'bo', 'linewidth', 1.5);
plot(0.945583025550488, 0.894127258067609, 'rx', 'linewidth', 1.4);
title('Final convergence');
xlabel('x_1'); ylabel('x_2');
xlim([0.6 1.3]); ylim([0 4.5]);
legend('Barrier','Convergence path','Iterate','Constrained minimum','Location','SouthEast')

%zoom
figure
plot(xaux,fxp, 'r', 'linewidth', 1); grid on; hold on
plot(xrvec(1,:),xrvec(2,:), 'g','linewidth', 1.5); 
plot(xrvec(1,:),xrvec(2,:), 'bo', 'linewidth', 1.5);
plot(0.945583025550488, 0.894127258067609, 'rx', 'linewidth', 1.4);
title('Final convergence (zoom)');
xlabel('x_1'); ylabel('x_2');
xlim([xP(1)-0.001 xP(1)+0.001]); ylim([xP(2)-0.002 xP(2)+0.002]);
legend('Barrier','Convergence path','Iterate','Constrained minimum','Location','SouthEast')

%%

if exist('fmincon') ~= 0                                  % MATLAB Optimization Toolbox
    options = optimoptions(@fmincon,'Algorithm','sqp');
    [xmin1,fval, exitflag,output] = fmincon(fx,x0,[],[],[],[],[],[],... 
       @confuneqBarreira,options);
else                                                      % Octave: built-in sqp
    [xmin1, fval] = sqp(x0, fx, @(x) x(1)^2-x(2));
end
fprintf('Reference (SQP): x = [%.4f; %.4f]\n', xmin1);



