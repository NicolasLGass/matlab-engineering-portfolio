%% Metodo Barreira
function [xi,xrvec,i,vecval, vecalfa] = MBarreira(fx,g,B,dB,x0,ci,gamma,x1,x2,nim,errp,itmaxDico,results)
xi= x0;

phi= @(x,ci) fx(x) - ci*B(x);                   %Função penalizada
dphi= @(x,ci) g(x) - ci*dB(x);                  %Derivada da função phi
    
i= 0;
normi= 10;                                      %Norma ini diff xi e xi+1

xrvec= xi;                                      %Vetor valores de xi

while normi>errp && i<=nim
    
    [xt, vecval, vecalfa]= MQsNewtDicoB(phi,dphi,B,xi,ci,x1,x2,nim,errp,itmaxDico,results);
%     [xt, vecval, vecalfa] = BuscaMulti_QuaseNewton(phi,dphi,B,ci,xi,tol,1500,'DFP',1);
%     phixi= phi(xi,ci)
%     phixt= phi(xt,ci)
%     xt
%     ci
    
    
    normi= norm(xt-xi,2);                       %Calc variação do pt xi
    if normi>errp
        ci= ci/gamma;                               %Recalculando gamma
    end
    
%     normi
%     pause
    
    xi= xt;                                     
    i= i+1;
    
    xrvec= [xrvec, xi];
    
    if results== 1
        fprintf('Iterações: %d  ci: %.8f \n',i,ci);
        fprintf('xi: %.4f ; %.4f ||xi+1 - xi||: %.6f \n',xi(1),xi(2),normi);
        
%         figure(2)
%         plot(xrvec(1,:),xrvec(2,:), 'g','linewidth', 1.5); grid on; hold on
%         plot(xrvec(1,:),xrvec(2,:), 'bo', 'linewidth', 1.5);
% %         plot(0.25,0.75, 'rx', 'linewidth', 1.4);
%         title('Convergência');
%         xlabel('Valor de x1'); ylabel('Valor de x2');
%         xlim([0 0.26]); ylim([0 0.78]);
%         legend('Convergência','Valor da iteração','Mínimo absoluto','Location','SouthEast')
        pause
    end
    
end
end
