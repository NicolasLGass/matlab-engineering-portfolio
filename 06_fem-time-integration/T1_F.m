%% Elementos Finitos 3 - Trabalho 1 - Elemento triangular
%                                     Matriz de rigidez e massa
clear all; close all; clc;
if exist('OCTAVE_VERSION', 'builtin'), pkg load symbolic; end   % Octave: symbolic package

%% Parametros iniciais

syms A x y x1 x2 x3 y1 y2 y3;

N1= (1/(2*A))* ((x2*y3-x3*y2)+(y2-y3)*x+(x3-x2)*y);
N2= (1/(2*A))* ((x3*y1-x1*y3)+(y3-y1)*x+(x1-x3)*y);
N3= (1/(2*A))* ((x1*y2-x2*y1)+(y1-y2)*x+(x2-x1)*y);

N= [N1 0 N2 0 N3 0;...
    0 N1 0 N2 0 N3];

dNx= diff(N,x);
dNy= diff(N,y);

%% Matriz de rigidez

% Derivadas de N sao constantes no elemento: a integral de area vira produto por A
Ke1= (dNx.'*dNx)*A;
Ke2= (dNy.'*dNy)*A;

disp('Matriz de rigidez:')
Ke= Ke1+Ke2                                     % Matriz de rigidez

%% Matriz de massa

syms rho te xi eta;

% Integral de area no triangulo: mapeamento para o triangulo de referencia
% x = x1+(x2-x1)*xi+(x3-x1)*eta, y = y1+(y2-y1)*xi+(y3-y1)*eta, dA = 2A dxi deta
Nref= subs(N, [x y], [x1+(x2-x1)*xi+(x3-x1)*eta, y1+(y2-y1)*xi+(y3-y1)*eta]);

disp('Matriz de massa ("A" simbólico):')
me= rho*te*int(int(Nref.'*Nref*2*A, xi, 0, 1-eta), eta, 0, 1)   % Matriz de massa ('A' simbólico)

At= [1 x1 y1;...
     1 x2 y2;...
     1 x3 y3];
An= (1/2)*det(At);

disp('Matriz de massa ("A" substituido):')
mef= simplify(subs(me,A,An))                    % Matriz de massa final (= rho*te*An/12*[2 1 1 ...])

%% Matriz de massa concentrada (metodo soma de linhas)

m1= mef(1,1)+mef(1,2)+mef(1,3)+mef(1,4)+mef(1,5)+mef(1,6);
m2= mef(2,1)+mef(2,2)+mef(2,3)+mef(2,4)+mef(2,5)+mef(2,6);
m3= mef(3,1)+mef(3,2)+mef(3,3)+mef(3,4)+mef(3,5)+mef(3,6);
m4= mef(4,1)+mef(4,2)+mef(4,3)+mef(4,4)+mef(4,5)+mef(4,6);
m5= mef(5,1)+mef(5,2)+mef(5,3)+mef(5,4)+mef(5,5)+mef(5,6);
m6= mef(6,1)+mef(6,2)+mef(6,3)+mef(6,4)+mef(6,5)+mef(6,6);

dm=[m1 m2 m3 m4 m5 m6];

disp('Matriz de massa concentrada (metodo soma de linhas, não normalizada):')
mlumped= diag(dm)                               % Matriz de massa concentrada

%% Matriz de massa concentrada HRZ

sumtot= 0;                                      % Inicializanto variaveis de somatorio
sumdg= 0;
sumndg= 0;
contv= 1;                                       % cont variaveis fora da diagonal
contvd= 1;                                      % cont variaveis na diagonal

for i=1:6
    for j=1:6
        sumtot= sumtot+me(i,j);                 % Soma da massa total
        if i~=j
            sumndg= sumndg+me(i,j);             % Soma da massa fora da diagonal
            fdiag(contv)= me(i,j);
            contv= contv+1;
        else
            sumdg= sumdg+me(i,j);               % Soma da massa da diagonal
            dg(contvd)= me(i,j);
            contvd= contvd+1;
        end
    end
end

totpndg= sumtot/sumndg;                         % Razão entre massa total e massa fora da diagonal

dgf= dg/totpndg;                                % Dividindo a diagonal pela razão

disp('Matriz de massa concentrada (HRZ):')
MHRZ= diag(dgf)                                 % Matriz de massa concentrada HRZ
