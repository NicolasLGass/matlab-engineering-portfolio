%% Método do gradiente com busca dicotômica
% f- função avaliada
% x0- ponto inicial
% x1- ponto inicial intervalo de analise (Dicotômica)
% x2- ponto final intervalo de analise (Dicotômica)
% nim- número max de iterações
% errp- erro mínimo aceito
% itmaxDico- Iterações máximas do método dicotômico
% results- mostrar resultados dos métodos (0 ou 1)

function [val, vecalfa] = MQsNewtDicoP(phi,dphi,xi,DadosIni)
warning off;

%Estrutura de dados inicias
ci= DadosIni.ci;  gci= DadosIni.gci;
h= DadosIni.h;  hj= DadosIni.hj;  ghj= DadosIni.ghj;      %Par. Lagrangeano
nim= DadosIni.nim;  errp= DadosIni.errp;
results= DadosIni.results;

lambj= DadosIni.lamb0;  nui= DadosIni.nu0;
ri= DadosIni.ri0;  rj= DadosIni.rj0;

Nx0= size(xi,1);                            %Tamanho x0
I= eye(Nx0);                                %Matriz identidade
Di= I;                                      %Di= Matriz identidade

%Definindo método DFP
Daux= @(Dn,yn,sn) Dn+(sn*sn.')/(sn.'*yn)-(Dn*(yn*yn.')*Dn)/(yn.'*Dn*yn);

val= xi;                                    %Valor inicial
vecval= [val];                               %Vetor de valores
vecalfa= [0];

for k=1:1:nim

    di= -Di*dphi(val,nui,lambj,ri,rj,hj,ghj);                 

    alfa= MDicoP(phi,val,di,DadosIni,errp/100,errp,0); %Dicotómica
    
    valalp= val+alfa*di;                  %Calculado Valor de alpha(x(i+1))

    %Calculo y e s do método DFP
    s= valalp-val;
    y= dphi(valalp,nui,lambj,ri,rj,hj,ghj)-dphi(val,nui,lambj,ri,rj,hj,ghj);
    
    %Calculo das normas
	dx1= norm(s,2);
    dx2= norm(phi(valalp,nui,lambj,ri,rj,h,hj)-phi(val,nui,lambj,ri,rj,h,hj),2);
        
    if s.'*y>0                          %teste se a matriz é positiva
        Di= Daux(Di,y,s);
    else
        Di= I;
    end
    
    val= valalp;                        %Atualiza o valor de xi
    vecval= [vecval, val];              %Atualiza vetor de valores
    
    if errp >= dx1 && errp >= dx2
        break;
    end
end %for k

%plot opcional de resultados
if results==1
    figure(1)
    plot(vecval(1,:),vecval(2,:), 'bo', 'linewidth', 1.2); hold on; grid on;
    plot(vecval(1,:),vecval(2,:), 'g', 'linewidth', 1);
    plot(0.25,0.75, 'rx', 'linewidth', 1.4);
    title('Iterates until convergence - augmented Lagrangian method');
    xlabel('x_1'); ylabel('x_2');
    xlim([0 0.42]); ylim([0 0.78]);
    legend('Iterate','Convergence path','Global minimum','Location','SouthEast')    
end %if results
end