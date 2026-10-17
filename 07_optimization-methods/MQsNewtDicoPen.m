%% Metodo Quase-Newton (DFP) com busca dicotomica para funcao penalizada phi(x,ci)
% Subproblema irrestrito do Metodo de Penalidade (MPenalidade.m).
% Variaveis de entrada:
% phi, dphi - Funcao penalizada e seu gradiente, phi(x,ci) e dphi(x,ci)
% x0        - Ponto inicial
% ci        - Parametro de penalidade atual
% x1, x2    - Intervalo do passo alfa para a busca dicotomica
% nim       - Numero maximo de iteracoes
% errp      - Tolerancia (variacao de x e de phi)
% itmaxDico - Numero maximo de iteracoes da busca dicotomica
% results   - (nao utilizado; mantido pela assinatura)
% Saidas:
% val       - Minimo de phi(.,ci)
% vecval    - Historico dos pontos

function [val, vecval] = MQsNewtDicoPen(phi,dphi,x0,ci,x1,x2,nim,errp,itmaxDico,results) %#ok<INUSD>

Nx0= size(x0,1);
I= eye(Nx0);
Di= I;                                      %Di inicial = identidade
Daux= @(Dn,yn,sn) Dn+(sn*sn.')/(sn.'*yn)-(Dn*(yn*yn.')*Dn)/(yn.'*Dn*yn);    % Atualizacao DFP

val= x0;
vecval= x0;
TolD= 0.001*(x2-x1);

for k=1:nim
    di= -Di*dphi(val,ci);
    alfa= MDicoPen(phi,ci,val,x1,x2,di,TolD/10,TolD,itmaxDico,0);
    valalp= val+alfa*di;

    s= valalp-val;
    y= dphi(valalp,ci)-dphi(val,ci);
    dx1= norm(s,2);
    dx2= abs(phi(valalp,ci)-phi(val,ci));

    if s.'*y>0                              %Mantem Di positiva definida
        Di= Daux(Di,y,s);
    else
        Di= I;
    end
    val= valalp;
    vecval= [vecval, val]; %#ok<AGROW>

    if errp >= dx1 && errp >= dx2
        break;
    end
end
end
