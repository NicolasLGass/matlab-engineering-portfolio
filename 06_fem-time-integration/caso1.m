% Dados iniciais primeiro caso, equilatero

% Posição dos nos
%    n   x2  y3
nos=[1   0   0;
     2   3   4;
     3   6   0];

% Nos de conexao dos elementos
%   n   i2  f3
el=[1   1   2;
    2   2   3];

% Relação (A*E)/L
eal= 1;

% Forcas aplicadas
%  n  noh XouY? F
F=[1   2   1   1];

% Restricoes:
%    noh XouY valor
Restr=[1   1   0
       1   2   0
       3   1   0
       3   2   0];
   

