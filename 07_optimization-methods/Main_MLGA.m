% Main Metodo da Penalidade
clear all; close all; clc;

% Função P1

fx= @(x) 1/2*(x(1)^2 + (x(2)^2)/3);         %Função inicial
g= @(x) [x(1); x(2)/3];                     %Grandiente da função

ci= @(x) (x(1)+x(2)-1);                     %Penalização
gci= @(x) [1; 1];                           %Gradiente da penalização

x0= [0; 0];                                 %x inicial
xf= [0.25; 0.75];                           %x ótimo final

h= 0;
hj= 0;
ghj= [0; 0];

lamb0= 0;                                   %Multiplicador de Lagrange Lambda
nu0= 0;                                     %Multiplicador de Lagrange Nu
ri0= 1;                                     %Par. de penalidade iniciais
rj0= 1;

errp= 0.00001;                              %Erro de tolerância
nim= 100;                                   %Número máx de iterações
itmax= nim;                                 %N max de iterações dicotômica
results= 1;                                 %Mostrar resultados, 1=sim

x1=0;                                       %Intervalo dicotômica 
x2=1;

%Estrutura de dados inicias
DadosIni.fx= fx;  DadosIni.g= g;                         %Funções
DadosIni.ci= ci;  DadosIni.gci= gci;
DadosIni.x0= x0;
DadosIni.h= h;  DadosIni.hj= hj;  DadosIni.ghj= ghj;      %Par. Lagrangeano
DadosIni.lamb0= lamb0;  DadosIni.nu0= nu0;               
DadosIni.ri0= ri0;  DadosIni.rj0= rj0;                   %Par. Penalidade
DadosIni.x1= x1;  DadosIni.x2= x2;                       %Int dicotomica
DadosIni.nim= nim;  DadosIni.errp= errp;  DadosIni.itmax= itmax;
DadosIni.results= results;

[xP,xrvec,i]= MLGAum(DadosIni);

figure
plot(xrvec(1,:), 'linewidth', 1.5); grid on; grid minor;
title('Convergence of x_1');
xlabel('Iteration'); ylabel('x_1');
ylim([0 0.26])

figure
plot(xrvec(2,:), 'linewidth', 1.5); grid on; grid minor;
title('Convergence of x_2');
xlabel('Iteration'); ylabel('x_2');

%% fmincon
if exist('fmincon') ~= 0                                  % MATLAB Optimization Toolbox
    options = optimoptions(@fmincon,'Algorithm','sqp');
    [xmin1,fval, exitflag,output] = fmincon(fx,x0,[],[],[],[],[],[],... 
       @cfeqLGA,options);
else                                                      % Octave: built-in sqp
    [xmin1, fval] = sqp(x0, fx, @(x) x(1)+x(2)-1);
end
fprintf('Reference (SQP): x = [%.4f; %.4f]\n', xmin1);
