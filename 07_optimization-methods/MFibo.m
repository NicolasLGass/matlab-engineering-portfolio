%% Método de Fibonacci
% Método que utiliza a série de Fibonacci para gerar relações de valores
% entre os valores dos limites de análise da função com os valores da série
% de Fibonacci, para encontrar o mínimo local da função

function MFibo(f,x1,x2,itmax)
format short;

%Sério de Fibonacci
fibv= ones(1,itmax);                          %Sério de Fibonacci
for i=3:itmax+1
    fibv(i)= fibv(i-1)+fibv(i-2);
end

fprintf('Resultados \n');
fprintf('Iteração  |   x1   |   x2   |   lmb   |   xm   |  mi  |  f(lmb)  |  f(xm)  |  f(mi) \n');
% fprintf('   %d       %.4f   %.4f   %.4f   %.4f   %.4f   %.4f     %.4f    %.4f \n', k, x1, x2, lmbk, xoti, mik, flmb, foti, fmi);
for k=1:itmax
    
%     lmbk= x2+((fibv(itmax-k-1))./(fibv(itmax-k+1))).*(x2-x1);
    mik= x1+((fibv(itmax-k+1))./(fibv(itmax-k+2))).*(x2-x1);
    %or
    lmbk= x1+x2-mik;
    
    flmb= f(lmbk);
    fmi= f(mik);
    xoti= (x1+x2)/2;                            %x ótimo
    foti= f(xoti);                              %Valor ótimo na função
    
    fprintf('  %d       %.4f   %.4f   %.4f   %.4f   %.4f   %.4f     %.4f    %.4f \n', k, x1, x2, lmbk, xoti, mik, flmb, foti, fmi);
    
    clf                     %Fecha as figuras anteriores
    figure
    fplot(f,[x1 x2], 'k', 'linewidth', 2); hold on ;grid on;             %Plot da função
    plot(xoti, foti, 'gx', 'linewidth', 2);           %Plot do Pt Médio
    plot(lmbk, flmb, 'rx', 'linewidth', 2);           %Plot do Pt Médio
    plot(mik, fmi, 'bx', 'linewidth', 2);           %Plot do Pt Médio
    legend('Function', sprintf('Midpoint %d', k), 'Point \lambda', 'Point \mu');
    title(sprintf('Iteration %d', k));
    pause 
    
    if flmb<fmi
        x2= mik;
    elseif flmb>fmi
        x1= lmbk;
    elseif flmb==fmi
        if min(abs(mik),abs(x1))==abs(x1)
            x2= lmbk;
        else
            x1= mik;
        end
    end
    
end %For

end %Function
