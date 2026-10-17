%% Método do gradiente com busca dicotômica
% f- função avaliada
% x0- ponto inicial
% x1- ponto inicial (Dicotômica)
% x2- ponto final (Dicotômica)
% nim- número max de iterações
% errp- erro mínimo aceito
% itmaxDico- Iterações máximas do método dicotômico
% results- mostrar resultados dos métodos (0 ou 1)

function [val, vecalfa] = MQsNewtDico(f,x0,x1,x2,nim,errp,itmaxDico,results)
warning off;

Nx0= size(x0,1);
I= eye(Nx0);
Di= I;

%Definindo método DFP
Daux= @(Dn,yn,sn) Dn+ (sn*sn.')/(sn.'*yn)- (Dn*(yn*yn.')*Dn)/(yn.'*Dn*yn);

val= x0;                                    %Valor inicial
vecval= [x0];                               %Vetor de valores
vecalfa=[0];

tic
for k=1:1:nim
    [~, gi]= f([val]);
    di= -Di*gi;
    
    alfa= MDicoF(f,val,x1,x2,di,errp/100,errp,itmaxDico,0); %Dicotómica
    vecalfa= [vecalfa, alfa];               %Vetor de Alfa(teste de conv)
    
    valalp= val+alfa*di;

    %Calculo y e s do método DFP
    [~, gval]= f([val]);                    %Gradiente de val
    [~, galp]= f([valalp]);                 %Gradiente de valalp
    y= galp-gval;
    s= valalp-val;
    
    if s.'*y>0                              %teste se a matriz é positiva
        Di= Daux(Di,y,s);
    else
        Di= I;
    end
    
    val= valalp;                        %Atualiza o valor de xi
    vecval= [vecval, val];              %Atualiza vetor de valores
    
    [~, gt]= f([val]);      % VER SE CONDICIONAL DE PARADA DEVE IR NO INICIO
	dx= norm(gt,2);
    if errp > dx
        fprintf('Convergência satisfeita... \n');
        fprintf('Iterações: %d \n',k);
        break;
    end
end %for k
toc
k;

%plot opcional de resultados
if results==1
    
    figure(1)
    plot(vecval(1,:),vecval(2,:), 'bo', 'linewidth', 1.2); hold on; grid on;
    plot(vecval(1,:),vecval(2,:), 'g', 'linewidth', 1);
    title('Iterates until convergence'); xlabel('x_1'); ylabel('x_2');
    incrp2=0.3;
    xlim([min(vecval(1,:))-incrp2 max(vecval(1,:))+incrp2])
    ylim([min(vecval(2,:))-incrp2 max(vecval(2,:))+incrp2])
%     datatip(-1.2,1);

    if f(x0)<10
        faux= @(x1,x2) 10*(x2-x1^2)^2+(1-x1)^2;     %Função auxiliar plot
    else
        faux= @(x1,x2) 200*(x2-x1^2)^2+(1-x1)^2;    %Função auxiliar plot
    end
    
    for i=1:length(vecval)
        vecx(i)= vecval(1,i);
        vecy(i)= vecval(2,i);
        vecz(i)= f([vecval(1,i),vecval(2,i)]);
    end
    
    % Calculando os máximos dos eixos do plot
    incrmt= 2;
    pxmin= min(vecx)-incrmt;
    pxmax= max(vecx)+incrmt;
    pymin= min(vecy)-incrmt;
    pymax= max(vecy)+incrmt;
    pzmin= min(vecz)-incrmt;
    pzmax= max(vecz)+incrmt;
    
    figure(4)
    [XXs, YYs] = meshgrid(linspace(-5, 5, 121));        % fsurf replacement (MATLAB/Octave)
    surf(XXs, YYs, arrayfun(faux, XXs, YYs), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on; grid on;
    plot3(vecval(1,1:end),vecval(2,1:end), vecz,'g', 'linewidth', 2);
    plot3(vecval(1,1:end),vecval(2,1:end), vecz,'bo', 'linewidth', 1.5);
    zlim([pzmin pzmax]); ylim([pymin pymax]); xlim([pxmin pxmax]);
    view(40,35);
%     datatip(faux,vecval(1,end),vecval(2,end),vecz(end));
end %if results
end