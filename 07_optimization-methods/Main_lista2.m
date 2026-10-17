%% Main Lista 2_2
% Lista E2 4.1.1; 4.1.2; 4.2; 4.3; 
% 4.1.1; 4.1.2 -> Implementar os métodos do gradiente, Newton
% 4.2 -> Verificar seus resultados através das funções fminsearch e
% fminunc do aplicativo MATLAB.
% 4.3 -> Realizar manualmente as 3 primeiras iterações de cada método
% para a função f1(x1, x2).
clear all; close all; clc;

%% Método do gradiente

%% Função f1 Método Gradiente
func= @func1;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x01= [-1.2; 1];                             %Valor inicial
x1=0; x2= 0.01;                             %Limites de Alfa
errp= 0.001;                                %Dif de passo de convergência
itmaxDico= 100;                              %Int max M Dicotomico
results= 1;                                 %Resultados M Dicotomico

[val] = MGradDicoF(func,x01,x1,x2,nim,errp,itmaxDico,results);

%Comparação
[veriff11,~,~,info11]= fminsearch(func,[-1.2; 1]);
[veriff12,~,~,info12]= fminunc(func,[-1.2; 1]);
it11= info11.iterations;
it12= info12.iterations;

%% Função f2 Método Gradiente
clear all; close all;

func= @func2;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x02= [-1.2; 1];                             %Valor inicial
x1=0; x2= 0.0015;                            %Limites de Alfa
errp= 0.02;                                  %Dif de passo de convergência
itmaxDico= 100;                              %Int max M Dicotomico
results= 1;                                 %Resultados M Dicotomico

[val] = MGradDicoF(func,x02,x1,x2,nim,errp,itmaxDico,results);

%Comparação
[veriff21,~,~,info21]= fminsearch(func,[-1.2; 1]);
[veriff22,~,~,info22]= fminunc(func,[-1.2; 1]);
it21= info21.iterations;
it22= info22.iterations;

%% Função f3 Método Gradiente
clear all; close all; 

func= @func3;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x03= [0; 0; 0; 0];                                %Valor inicial
x1=0; x2= 10;                               %Limites de Alfa
errp= 0.01;                                 %Dif de passo de convergência
itmaxDico= 20;                              %Int max M Dicotomico
results= 0;                                 %Resultados M Dicotomico

[val] = MGradDicoF(func,x03,x1,x2,nim,errp,itmaxDico,results);

%% Função f4 Método Gradiente
clear all; close all; 

func= @func4;                               %Mudando chamada da função
nim= 20000;                                 %Número máximo de iterações

x03= [0; 0; 0; 0];                                %Valor inicial
x1=0; x2= 1;                            %Limites de Alfa
errp= 0.01;                                 %Dif de passo de convergência
itmaxDico= 20;                              %Int max M Dicotomico
results= 0;                                 %Resultados M Dicotomico

[val] = MGradDicoF(func,x03,x1,x2,nim,errp,itmaxDico,results);

%% Método do Newton

%% Função f1 Método Newton
clear all; close all;

func= @func1;                               %Mudando chamada da função
nim= 20000;                                  %Número máximo de iterações

x01= [-1.2; 1];                             %Valor inicial
x1=0; x2= 0.1;                             %Limites de Alfa
errp= 0.001;                                %Dif de passo de convergência
itmaxDico= 100;                              %Int max M Dicotomico
results= 1;                                 %Resultados M Dicotomico

[val] = MNewtDicoF(func,x01,x1,x2,nim,errp,itmaxDico,results);

%Comparação
[veriff31,~,~,info31]= fminsearch(func,[-1.2; 1]);
[veriff32,~,~,info32]= fminunc(func,[-1.2; 1]);
it31= info31.iterations;
it32= info32.iterations;

%% Função f2 Método Newton
clear all; close all;

func= @func2;                               %Mudando chamada da função
nim= 20000;                                  %Número máximo de iterações

x02= [-1.2; 1];                             %Valor inicial
x1=0; x2= 0.1;                             %Limites de Alfa
errp= 0.001;                                %Dif de passo de convergência
itmaxDico= 100;                              %Int max M Dicotomico
results= 1;                                 %Resultados M Dicotomico

[val] = MNewtDicoF(func,x02,x1,x2,nim,errp,itmaxDico,results);

%Comparação
[veriff41,~,~,info41]= fminsearch(func,[-1.2; 1]);
[veriff42,~,~,info42]= fminunc(func,[-1.2; 1]);
it41= info41.iterations;
it42= info42.iterations;

%% Função f3 Método Newton
clear all; close all;

func= @func3;                               %Mudando chamada da função
nim= 20000;                                  %Número máximo de iterações

x03= [0; 0; 0; 0];                                %Valor inicial
x1=0; x2= 10;                            %Limites de Alfa
errp= 0.01;                                 %Dif de passo de convergência
itmaxDico= 20;                              %Int max M Dicotomico
results= 0;                                 %Resultados M Dicotomico

[val] = MNewtDicoF(func,x03,x1,x2,nim,errp,itmaxDico,results);

%% Função f4 Método Newton
clear all; close all;

func= @func4;                               %Mudando chamada da função
nim= 20000;                                  %Número máximo de iterações

x04= [0; 0; 0; 0];                                %Valor inicial
x1=0; x2= 1;                            %Limites de Alfa
errp= 0.01;                                 %Dif de passo de convergência
itmaxDico= 20;                              %Int max M Dicotomico
results= 0;                                 %Resultados M Dicotomico

[val] = MNewtDicoF(func,x04,x1,x2,nim,errp,itmaxDico,results);

