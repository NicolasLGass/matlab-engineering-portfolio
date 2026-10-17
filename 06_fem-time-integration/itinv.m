%% Elementos Finitos 3 - Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
% ---------------------------------------------
% --- Iteracao inversa (Inverse Power Method) ---
% ---------------------------------------------
% Resolve K*phi = lambda*M*phi para o menor autovalor ainda nao encontrado.
% Os modos ja obtidos (Phi) sao removidos a cada iteracao por
% ortogonalizacao de Gram-Schmidt em relacao a M (deflacao).
%
% Variaveis de entrada:
% K - Matriz de rigidez (nxn)
% M - Matriz de massa (nxn)
% xini - Vetor inicial (nx1)
% intc - Numero maximo de iteracoes
% tol - Tolerancia relativa no autovalor
% Phi - (opcional) Modos ja obtidos, M-ortonormais (nxm)
%
% Saidas:
% phi - Modo normalizado pela massa (phi'*M*phi = 1)
% lambda - Autovalor (omega^2)
% nit - Numero de iteracoes realizadas

function [phi, lambda, nit]= itinv(K, M, xini, intc, tol, Phi)

[mk,nk]= size(K);
[mm,nm]= size(M);
if mk~=nk || mm~=nm
    error('Matriz K ou M nao sao quadradas.');
end
if nargin < 6
    Phi= zeros(mk,0);
end

x0= xini;
lambda= 0;

for nit= 1:intc
    lambdaold= lambda;                  % salva autovalor antigo

    x0= x0 - Phi*(Phi'*M*x0);           % remove componentes dos modos ja obtidos
    xib= K\(M*x0);                      % x barra
    xib= xib - Phi*(Phi'*M*xib);        % Gram-Schmidt (M-ortogonalizacao)

    x0= xib/sqrt(xib'*M*xib);           % normalizacao pela massa
    lambda= x0'*K*x0;                   % quociente de Rayleigh (x0'*M*x0 = 1)

    err= abs(lambda-lambdaold)/abs(lambda);   % Erro de iteracao
    if err <= tol
        break;
    end
end

if err > tol
    warning('itinv: nao convergiu em %d iteracoes (erro = %.2e).', intc, err);
end

phi= x0;
end
