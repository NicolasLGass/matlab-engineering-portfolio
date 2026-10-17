%%  Método da Busca Dicotômica pela primeira derivada
% A função divide o intervalo da função avaliada, substituindo os valores
% do intervalo, para convergir e englobar o mínimo valor da função.
% func- função avaliada
% x0- ponto inicial
% x1- ponto inicial (Dicotômica)
% x2- ponto final (Dicotômica)
% di- vetor de direção
% erro- erro de tolerância
% intrvmin- intervalo mínimo aceito
% itmax- Iterações máximas
% results- mostrar resultados dos métodos (0 ou 1)

%                   (phi,val,di,DadosIni,errp/100,errp,0)
function xm = MDicoP(func,x0,di,DadosIni,erro,intrvmin,results)
warning off

%Estrutura de dados inicias
ci= DadosIni.ci;  hj= DadosIni.hj; h= DadosIni.h;     %Par. Lagrangeano
lambj= DadosIni.lamb0;  nui= DadosIni.nu0;
ri= DadosIni.ri0;  rj= DadosIni.rj0;
x1= DadosIni.x1;  x2= DadosIni.x2;
itmax= DadosIni.itmax;

k= 0;                                %Inicialização contador de iterações
intrv= (x2-x1);                     %Intervalo inicial

if results == 1
    fprintf('Resultados \n');
    fprintf('Iteração  |   x1   |   x2   |   lmb   |   xm   |  mi  |  f(lmb)  |  f(xm)  |  f(mi)   |   Intervalo  \n');
end

while intrv>intrvmin && k<=itmax
    
    xm= (x2+x1)/2;                      %Ponto Médio do intervalo
    lmb= xm-erro;                       %Ponto intermediário menor
    mi= xm+erro;                        %Ponto intermediário maior
    xm2= xm*di;                         %Ponto Médio do intervalo
    lmb2= x0+ lmb*di;
    mi2= x0+ mi*di;
    
    fval= func(xm2,nui,lambj,ri,rj,h,hj);          %Valor de f em xm
    flmb= func(lmb2,nui,lambj,ri,rj,h,hj);         %Valor de f em lambda
    fmi= func(mi2,nui,lambj,ri,rj,h,hj);           %Valor de f em mi
    
    k= k+1;
    
    if results == 1
        fprintf('  %d       %.4f   %.4f   %.4f   %.4f   %.4f   %.4f     %.4f    %.4f      %.4f \n',k, x1, x2, lmb, xm, mi, flmb, fval, fmi, intrv);
    end
    
    if flmb>fmi                      %Se for decrescente para a esquerda
        x1= lmb;                     %Troca o valor do lim superior
    else
        x2= mi;                      %Troca o valor do lim inferior
    end
    
    intrv= (x2-x1);                  %Intervalo
    
end %While
xm= (x2+x1)/2;                       %Ponto Médio do intervalo
end %Function
