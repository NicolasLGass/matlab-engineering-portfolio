%% Elementos Finitos 3 - Autovalores e modos pela iteracao inversa com deflacao
% Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
clear variables; close all; clc;

%% Sistema de 4 GDL

K= [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 1];    % Matriz de rigidez
M= diag([1 2 1 1]);                                 % Matriz de massa
intmax= 1000;                                       % Numero maximo de iteracoes
tol= 1e-12;                                         % Tolerancia

%% Iteracao inversa: um modo por vez, deflacionando os ja obtidos

n= length(M);
Phi= zeros(n,n);
lambda= zeros(n,1);
nit= zeros(n,1);
xini= ones(n,1);                                    % Vetor inicial

for i= 1:n
    [Phi(:,i), lambda(i), nit(i)]= itinv(K, M, xini, intmax, tol, Phi(:,1:i-1));
end

%% Validacao com eig

lambdaEig= sort(eig(K, M));

fprintf('Modo   lambda (itinv)      lambda (eig)      iteracoes\n');
for i= 1:n
    fprintf('%3d   %16.12f   %16.12f   %6d\n', i, lambda(i), lambdaEig(i), nit(i));
end

fprintf('\nOrtogonalidade Phi''*M*Phi (deve ser a identidade):\n');
disp(round(Phi'*M*Phi * 1e10) / 1e10);       % round(x, n) is MATLAB-only
