%% Lista 1 - Apostila 3.7
% Seja f(x) = x2 + 2x e um intervalo inicial I0 = [?3, 5]. Realizar 
% mente três iterações dos seguintes métodos de busca unidimensional:
% 1. Busca dicotômica; 2. Busca Fibonacci; 3. Busca de Newton.
clear all; close all; clc;
if exist('OCTAVE_VERSION', 'builtin'), pkg load symbolic; end   % Octave: symbolic package
%%

x1= -3;
x2= 5;
x= linspace(x1,x2,100);

f= @(x) x.^2+2*x;
% f= @(x) x.^4;

erro= 0.001;
intrvmin= 0.01;                     %Intervalo Mínimo 
itmax= 10;                           %Máximo de iterações (Fibonacci)
erro2= 0.001;                         %Erro Método de Newton
xk= 4;                              %Valor inicial de xk (Método de Newton)

%% Dicotômica padrão
MDicotomica2(f,x1,x2,erro,intrvmin,itmax,1);

%% Dicotômica pela primeira derivada
MDicotomicaDeriv(f,x1,x2,erro,intrvmin)
%% Fibonacci
MFibo(f,x1,x2,itmax)
%% Newton
results=0;
[xNewt, fxNewt]= MNewton(f,x1,x2,xk,erro2,itmax,results);
