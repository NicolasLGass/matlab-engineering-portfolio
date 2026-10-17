%% Metodo Lagrangeano Aumentado
function [xi,xrvec,i] = MLGAum(DadosIni)

%Estrutura de dados inicias
fx= DadosIni.fx;  g= DadosIni.g;                         %Funções
ci= DadosIni.ci;  gci= DadosIni.gci;
h= DadosIni.h;  hj= DadosIni.hj;  ghj= DadosIni.ghj;      %Par. Lagrangeano
nim= DadosIni.nim;  errp= DadosIni.errp;
results= DadosIni.results;

xi= DadosIni.x0;
lambj= DadosIni.lamb0;  nui= DadosIni.nu0;
ri= DadosIni.ri0;  rj= DadosIni.rj0;

phi= @(x,nui,lambj,ri,rj,h,hj) fx(x) + nui*ci(x) + lambj*h + 1/2*ri*(ci(x).^2) + 1/2*rj*(hj.^2);                   %Função penalizada
dphi= @(x,nui,lambj,ri,rj,hj,ghj) g(x) + nui.*gci(x) + lambj.*ghj + ri.*ci(x).*gci(x) + rj.*hj.*ghj;                  %Derivada da função phi

i= 0;
normi= 1;                                      %Norma ini diff xi e xi+1

xrvec= xi;                                      %Vetor valores de xi

while normi>errp && i<=nim
    
    xt= MQsNewtDicoP(phi,dphi,xi,DadosIni);
    
    normi= norm(xt-xi,2);                       %Calc variação do pt xi
    
    if abs(ci(xt))>errp                                %Atualizando ci
        DadosIni.nu0= DadosIni.nu0+DadosIni.ri0*(ci(xt));    % atualiza com o r atual
        DadosIni.ri0= 2*DadosIni.ri0;                         % penalidade dobra a cada iteracao
    end
    if hj>0                                 %Atualizando lambj
        DadosIni.lamb0= DadosIni.lamb0+DadosIni.rj0*(hj);
        DadosIni.rj0= 2*DadosIni.rj0;
    end

    xi= xt;
    i= i+1;
    
    xrvec= [xrvec, xi];
    
    if results==1
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
