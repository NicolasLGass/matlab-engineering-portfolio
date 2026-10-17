%%  Método da Busca Dicotômica pela primeira derivada
% A função divide o intervalo da função avaliada, substituindo os valores
% do intervalo, para convergir e englobar o mínimo valor da função.

function xm = MDicotomica2(f,x1,x2,erro,intrvmin,itmax,results)
warning off

if erro>intrvmin
    fprintf('Cuidado, valor de intervalo mínimo menor que valor de erro admitido... \n');
end

intrv= abs(x1-x2);                  %Intervalo inicial
xm= (x2+x1)/2;                      %Ponto Médio do intervalo
lmb= ((x1+x2)/2)-erro;              %Ponto intermediário menor
mi= ((x1+x2)/2)+erro;               %Ponto intermediário maior

fval= f(xm);                   %Valor de f em xm
flmb= f(lmb);                  %Valor de f em lambda
fmi= f(mi);                    %Valor de f em mi

xmvec = xm;

k=1;                                %Inicialização contador de iterações
cntxm=0;                            %Contador Xm repetido

if results == 1
    fprintf('Resultados \n');
    fprintf('Iteração  |   x1   |   x2   |   lmb   |   xm   |  mi  |  f(lmb)  |  f(xm)  |  f(mi)   |   Intervalo  \n');
    fprintf('  %d       %.4f   %.4f   %.4f   %.4f   %.4f   %.4f     %.4f    %.4f      %.4f \n',k, x1, x2, lmb, xm, mi, flmb, fval, fmi, intrv);
end

if xm==0
    xmold=1;
elseif xm~=0
    xmold=0;
end

% while intrv>intrvmin
while abs(xmold-xm) > erro
    
    if results == 1
        clf                     %Fecha as figuras anteriores
        figure
        fplot(f,[x1 x2], 'k', 'linewidth', 2); hold on ;grid on;             %Plot da função
        plot(xm, fval, 'gx', 'linewidth', 2);           %Plot do Pt Médio
        plot(lmb, flmb, 'rx', 'linewidth', 2);           %Plot do Pt Médio
        plot(mi, fmi, 'bx', 'linewidth', 2);           %Plot do Pt Médio
        legend('Function', sprintf('Midpoint %d', k), 'Point \lambda', 'Point \mu');
        title(sprintf('Iteration %d', k));
        pause 
    end
    
    if flmb>fmi                      %Se df for decrescente para a esquerda
        x1= lmb;                     %Troca o valor do lim superior
    end
    if flmb<fmi                      %Se df for decrescente para a direita
        x2= mi;                      %Troca o valor do lim inferior
    end
%     if dfval==0      %Se df for decrescente para a direita
%         fprintf('Encerrado, ponto xm de valor %.4f é mínimo local \n', xm);
%         break
%     end    

    xmold=xm;
    
    intrv= abs(x1-x2);                  %Intervalo inicial
    xm= (x2+x1)/2;                      %Ponto Médio do intervalo
    lmb= (x1+x2)/(2)-erro;              %Ponto intermediário menor
    mi= (x1+x2)/(2)+erro;               %Ponto intermediário maior

    fval= f(xm);                   %Valor de f em xm
    flmb= f(lmb);                  %Valor de f em lambda
    fmi= f(mi);                    %Valor de f em mi
    
    xmvec= [xmvec, xm];
    
    k= k+1;
    
    if results == 1
        fprintf('  %d       %.4f   %.4f   %.4f   %.4f   %.4f   %.4f     %.4f    %.4f      %.4f \n',k, x1, x2, lmb, xm, mi, flmb, fval, fmi, intrv);
    end
    
    if k>=itmax
        fprintf('Valor máximo de iterações alcançado no método dicotómico.')
        fprintf('Valor de xm encontrado: %.4f \n', xm);
        break
    end
    
    if abs(xmold-xm) < erro
%         abs(xmold-xm)
%         fprintf('Valor de xm convergiu (<erro). \n')
%         fprintf('Valor de xm encontrado: %.4f \n', xm);
        break
    end
    
    if xmold==xm
        cntxm=cntxm+1;
        if cntxm>2
            fprintf('Valor de xm convergiu.')
            fprintf('Valor de xm encontrado: %.4f \n', xm);
            break
        end
    end
    
end %While

if intrv<intrvmin
    fprintf('Encerrado, intervalo menor que valor mínimo definido. \n');
    fprintf('Valor de xm encontrado: %.4f \n', xm);
end

end %Function