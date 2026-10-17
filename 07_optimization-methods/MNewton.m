%% Método de Newton

function [xk, fval]= MNewton(f,x1,x2,xk,erro2,itmax,results)
warning off
syms x

df= diff(f(x));                     %Derivada da função f
df2= diff(f(x),2);                  %Segunda derivada da função f
fval= f(xk);                %Valor de f em xm

dfxk= double(subs(df,x,xk));                  %Valor da derivada de f em xk
dfxk2= double(subs(df2,x,xk));                %Valor da segunda derivada de f em xk
rxk= -(dfxk)/dfxk2;                  %Relação de x-xk

xk1= xk+rxk;

k=1;

if results == 1
    fprintf('Resultados \n');
    fprintf('Iteração  |   xk  |  df(xk)  |  df2(xk)  |  xk+1 \n');
    fprintf('  %d         %.4f   %.4f   %.4f    %.4f    \n',k, xk, dfxk, dfxk2, xk1);

    clf                        %Fecha as figuras anteriores
    figure
    fplot(f,[x1 x2], 'k', 'linewidth', 2); hold on ;grid on;
    plot(xk, fval, 'rx', 'linewidth', 2);
    legend('Function', sprintf('x_k (iteration %d)', k));
    title(sprintf('Iteration %d', k));
    pause 
end

while abs(rxk)>erro2
    
    if k>=itmax
        break;                          %Condicional para máximo de iterações
    end
    
    xk=xk1;
    
    fval= f(xk);                %Valor de f em xm
    dfxk= double(subs(df,x,xk));                  %Valor da derivada de f em xk
    dfxk2= double(subs(df2,x,xk));                %Valor da segunda derivada de f em xk
    rxk= -(dfxk)/dfxk2;                 %Relação de x-xk
    xk1= xk+rxk;
    
    k= k+1;
    
    if results == 1

        fprintf('  %d       %.4f   %.4f   %.4f     %.4f    \n',k, xk, dfxk, dfxk2, xk1);

        clf                        %Fecha as figuras anteriores
        figure
        fplot(f,[x1 x2], 'k', 'linewidth', 2); hold on ;grid on;
        plot(xk, fval, 'rx', 'linewidth', 2);
        legend('Function', sprintf('x_k (iteration %d)', k));
        title(sprintf('Iteration %d', k));

        pause 
    end
    
end %While

end %Function