%% Metodo dicotomico para busca linear em funcao penalizada phi(x,ci)
% Usado pelos metodos de Penalidade (MQsNewtDicoPen) e Barreira (MQsNewtDicoB).
% Variaveis de entrada:
% func     - Funcao penalizada phi(x,ci)
% ci       - Parametro de penalidade/barreira atual
% x0       - Ponto atual
% x1, x2   - Intervalo inicial do passo alfa
% di       - Direcao de busca
% erro     - Distancia entre os pontos intermediarios (lmb, mi)
% intrvmin - Tamanho minimo do intervalo (criterio de parada)
% itmax    - Numero maximo de iteracoes
% results  - (nao utilizado; mantido pela assinatura)
% Saida:
% xm       - Passo alfa que minimiza phi(x0 + alfa*di, ci) em [x1, x2]

function xm = MDicoPen(func,ci,x0,x1,x2,di,erro,intrvmin,itmax,results) %#ok<INUSD>

k= 0;                                   %Inicialização contador de iterações
intrv= (x2-x1);                         %Intervalo inicial

while intrv>intrvmin && k<=itmax
    xm= (x2+x1)/2;                      %Ponto Médio do intervalo
    lmb= xm-erro;                       %Ponto intermediário menor
    mi= xm+erro;                        %Ponto intermediário maior

    flmb= func(x0+lmb*di, ci);          %Valor de phi em lambda
    fmi= func(x0+mi*di, ci);            %Valor de phi em mi
    k= k+1;

    if flmb>fmi                         %Se for decrescente para a direita
        x1= lmb;
    else
        x2= mi;
    end
    intrv= (x2-x1);
end

xm= (x2+x1)/2;                          %Ponto Médio do intervalo
end
