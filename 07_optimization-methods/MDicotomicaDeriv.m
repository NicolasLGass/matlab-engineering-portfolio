%%  Método da Busca Dicotômica Padrão
% A função divide o intervalo da função avaliada, substituindo o valor do
% ponto médio (xm) na derivada da função de entrada, assim indicando para
% o crescimento ou decrescimento da função para as direções avaliadas, com
% isso são substituídos os limites do intervalo para convergir e englobar o
% mínimo valor da função.

function MDicotomicaDeriv(f,x1,x2,~,intrvmin)
warning off
syms x

intrv= abs(x1-x2);                  %Intervalo inicial
xm= (x2+x1)/2;                      %Ponto Médio do intervalo
df= diff(f(x));                        %Derivada da função f
fval= f(xm);                   %Valor de f em xm
dfval= double(subs(df,x,xm));                 %Valor da derivada df em xm

k=1;

fprintf('Resultados \n');
fprintf('Iteração  |   x1   |   x2   |   xm   |  Intervalo  |  df(xm)  |  f(xm)  \n');
fprintf('  %d       %.4f   %.4f   %.4f      %.4f       %.4f     %.4f \n',k, x1, x2, xm, intrv, dfval, fval);

while intrv>intrvmin
    
    clf                     %Fecha as figuras anteriores
    figure
    fplot(f,[x1 x2], 'k', 'linewidth', 2); hold on ;grid on;
    plot(xm, fval, 'rx', 'linewidth', 2);
    legend('Function', sprintf('Midpoint %d', k));
    title(sprintf('Iteration %d', k));
    pause 
    
    if dfval>=0                      %Se df for decrescente para a esquerda
        x2= xm;                      %Troca o valor do lim superior
    end
    if dfval<0                       %Se df for decrescente para a direita
        x1= xm;                      %Troca o valor do lim inferior
    end
%     if dfval==0                      %Se df for decrescente para a direita
%         fprintf('Encerrado, ponto xm de valor %.4f é mínimo local \n', xm);
%         break
%     end
    
    intrv= abs(x1-x2);               %Recalculando Intervalo
    xm= (x2+x1)/2;                   %Recalculando Ponto Médio do intervalo
%     df= diff(f);                     %Derivada da função f
    fval= f(xm);                %Valor de f em xm
    dfval= double(subs(df,x,xm));              %Valor da derivada df em xm
    k= k+1;
    fprintf('  %d       %.4f   %.4f   %.4f      %.4f       %.4f     %.4f \n',k, x1, x2, xm, intrv, dfval, fval);
    
end %While

if intrv<intrvmin
    fprintf('Encerrado, intervalo menor que valor mínimo \n');
end

end %Function