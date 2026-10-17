%% Metodo Penalidade
function [xi,xrvec,i] = MPenalidade(fx,g,P,dP,x0,ci,gamma,x1,x2,nim,errp,itmaxDico,results)
xi= x0;

phi= @(x,ci) fx(x) + ci*P(x);
dphi= @(x,ci) g(x) + ci*dP(x);

i= 0;
normi= 10;

xrvec= xi;

while normi>errp && i<=nim
    
    xt= MQsNewtDicoPen(phi,dphi,xi,ci,x1,x2,nim,errp,itmaxDico,results);   % parte do ponto atual
    
    normi= norm(xt-xi,2);
    
    ci= ci*gamma;
    xi= xt;
    i= i+1;
    
    xrvec= [xrvec, xi];
    
    if results== 1
        fprintf('Iterações: %d  ci: %.4f \n',i,ci);
        fprintf('xi: %.4f ; %.4f ||xi+1 - xi||: %.6f \n',xi(1),xi(2),normi);
        
        figure(2)
        plot(xrvec(1,:),xrvec(2,:), 'g','linewidth', 1.5); grid on; hold on
        plot(xrvec(1,:),xrvec(2,:), 'bo', 'linewidth', 1.5);
        plot(0.25,0.75, 'rx', 'linewidth', 1.4);
        title('Convergence path');
        xlabel('x_1'); ylabel('x_2');
        xlim([0 0.26]); ylim([0 0.78]);
        legend('Convergence path','Iterate','Global minimum','Location','SouthEast')
        pause
    end
    
end
end
