%% Main Lista 3
% Lista E2 4.1.1; 4.1.2; 4.2; 4.3; 
% Exercício 4.4. Implementar o método da busca unidimensional da secante
% com o método Quase-Newton (DFP) e comparar os resultados com a implementa
% ção usando a busca dicotômica. Testar para as 4 funções anteriores.
clear all; close all; clc;

%% Função f1 Método Quase Newton
func= @func1;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x01= [-1.2; 1];                             %Valor inicial
x1=0; x2= 1;                                %Limites de Alfa
errp= 0.0001;                               %Dif de passo de convergência
itmaxDico= 1000;                            %Int max M Dicotomico
results= 1;                                 %Resultados M Dicotomico

[val, vecalfa] = MQsNewtDico(func,x01,x1,x2,nim,errp,itmaxDico,results);

%Comparação
% [veriff11,~,~,info11]= fminsearch(func,[-1.2; 1]);
% [veriff12,~,~,info12]= fminunc(func,[-1.2; 1]);
% it11= info11.iterations;
% it12= info12.iterations;

%% Função f2 Método Gradiente
clear all; close all;

func= @func2;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x02= [-1.2; 1];                             %Valor inicial
x1=0; x2= 1;                                %Limites de Alfa
errp= 0.0001;                                  %Dif de passo de convergência
itmaxDico= 100;                              %Int max M Dicotomico
results= 1;                                 %Resultados M Dicotomico

[val] = MQsNewtDico(func,x02,x1,x2,nim,errp,itmaxDico,results);

%Comparação
% [veriff21,~,~,info21]= fminsearch(func,[-1.2; 1]);
% [veriff22,~,~,info22]= fminunc(func,[-1.2; 1]);
% it21= info21.iterations;
% it22= info22.iterations;

%% Função f3 Método Gradiente
clear all; close all; 

func= @func3;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x03= [0; 0; 0; 0];                                %Valor inicial
x1=0; x2= 1;                               %Limites de Alfa
errp= 0.0001;                                 %Dif de passo de convergência
itmaxDico= 20;                              %Int max M Dicotomico
results= 0;                                 %Resultados M Dicotomico

[val] = MQsNewtDico(func,x03,x1,x2,nim,errp,itmaxDico,results);

%% Função f4 Método Gradiente
clear all; close all; 

func= @func4;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x04= [0; 0; 0; 0];                                %Valor inicial
x1=0; x2= 1;                            %Limites de Alfa
errp= 0.0001;                                 %Dif de passo de convergência
itmaxDico= 20;                              %Int max M Dicotomico
results= 0;                                 %Resultados M Dicotomico

[val] = MQsNewtDico(func,x04,x1,x2,nim,errp,itmaxDico,results);
