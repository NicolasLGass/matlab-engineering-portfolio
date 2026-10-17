%% Método do gradiente com busca dicotômica
% f- função avaliada
% x0- ponto inicial
% x1- ponto inicial intervalo de analise (Dicotômica)
% x2- ponto final intervalo de analise (Dicotômica)
% nim- número max de iterações
% errp- erro mínimo aceito
% itmaxDico- Iterações máximas do método dicotômico
% results- mostrar resultados dos métodos (0 ou 1)

function [val, vecval, vecalfa] = MQsNewtDicoB(f,dphi,B,x0,ci,x1,x2,nim,errp,itmaxDico,results)
warning off;

Nx0= size(x0,1);                            %Tamanho x0
I= eye(Nx0);                                %Matriz identidade
Di= I;                                      %Di= Matriz identidade

%Definindo método DFP
Daux= @(Dn,yn,sn) Dn+(sn*sn.')/(sn.'*yn)-(Dn*(yn*yn.')*Dn)/(yn.'*Dn*yn);

val= x0;                                    %Valor inicial
vecval= [x0];                               %Vetor de valores
vecalfa= [];
cont=0;

TolD = 0.001*(x2-x1);

for k=1:1:nim

    di= -Di*dphi(val,ci);
    
    alfa= MDicoPen(f,ci,val,x1,x2,di,TolD/10,TolD,itmaxDico,0); %Dicotómica
    
    valalp= val+alfa*di;                    %Calculado Valor de alpha(x(i+1))
    
    while B(valalp)>0                       %Verificando restrição
        alfa= 0.9*alfa;
        valalp= val+alfa*di;
%         B(valalp)
    end
	vecalfa= [vecalfa, alfa];
    
    %Calculo y e s do método DFP
    s= valalp-val;
    y= dphi(valalp,ci)-dphi(val,ci);

    %Calculo das normas
	dx1= norm(s,2);
    dx2= norm(f(valalp,ci)-f(val,ci),2);
%     pause
  
    if errp >= dx1 && errp >= dx2
        break;
    end
        
    if s.'*y>0                          %teste se a matriz é positiva
        Di= Daux(Di,y,s);
%         pause
    else
        Di= I;
    end
    
    val= valalp;                        %Atualiza o valor de xi
    vecval= [vecval, val];              %Atualiza vetor de valores
end %for k
k;

%plot opcional de resultados
if results==1
        
        xaux= 0.6:0.00001:2.05;
        for i=1:length(xaux)
            fxp(i)= (xaux(i).^2);
        end

        figure(1)
        plot(xaux,fxp, 'r', 'linewidth', 1); grid on; hold on
        plot(vecval(1,1),vecval(2,1), 'bo', 'linewidth', 1.5);
        plot(vecval(1,2:end-1),vecval(2,2:end-1), 'k*', 'linewidth', 1.2);
        plot(vecval(1,end),vecval(2,end), 'bo', 'linewidth', 1.5);
        plot(vecval(1,:),vecval(2,:), 'g', 'linewidth', 1);
        plot(0.945583025550488, 0.894127258067609, 'rx', 'linewidth', 1.6);
        title('Iterates until convergence');
        legend('Barrier','Final point','Intermediate points','Location','SouthEast')
        incrp2=0.1;

end %if results
end