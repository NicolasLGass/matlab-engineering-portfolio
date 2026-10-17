%% Elementos Finitos 3 - Analise flambagem de trelicas
% Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
clear all; close all; clc;
warning off;

%% Escolher caso da trelica - 1:Triangulo equilatero; 2:Triangulo retangulo
caso= 2;

if caso==1
    caso1;
else
    caso2;
end

% Plot da trelica em repouso
for i=1:size(el,1)
    x1=nos(el(i,2),2);
    x2=nos(el(i,3),2);
    y1=nos(el(i,2),3);
    y2=nos(el(i,3),3);
    figure(1)
    line([x1,x2],[y1,y2], 'Color', 'k', 'linestyle', '--', 'linewidth', 1.5); 
    hold on; grid on
end

%% Iniciando restricoes e parametros iniciais dos elementos

Nnos= size(nos,1);              % Numero de nos
GDLn= size(nos,2)-1;            % Numero de GDL por no
gdls= ones(Nnos,GDLn);          % Graus de liberdade de cada no

for i=1:size(Restr,1)
    Rno= Restr(i,1);            % No i de restricao
    Reixo= Restr(i,2);          % Eixo de restricao, x ou y
    Rval= Restr(i,3);           % Valore de restricao atribuido ao no
    gdls(Rno,Reixo)= Rval;      % Atribuindo o valore de restricao    
end

Ngdl=0;                         % Contador de GDL
for i=1:Nnos
    for j=1:GDLn
        if gdls(i,j)==1
            Ngdl= Ngdl+1;
            gdls(i,j)= Ngdl;    % Numerando os GDL
        end
    end
end

%% Matriz global

Keg= zeros(Ngdl,Ngdl);          % Matriz de rigidez global
Lvec= zeros(1,size(el,1));      % Comprimentos dos elementos (l)

for i=1:size(el,1)    
    no1= el(i,2);               % Primeiro no do elemento
    no2= el(i,3);               % Segundo no do elemento
	x1= nos(no1,2);             % Posicao em x do no 1 do elemento
    x2= nos(no2,2);             % Posicao em x do no 2 do elemento
    y1= nos(no1,3);             % Posicao em y do no 1 do elemento
    y2= nos(no2,3);             % Posicao em y do no 2 do elemento
    
    Lxi= x2-x1;                 % Comprimento do elemento no eixo x
    Lyi= y2-y1;                 % Comprimento do elemento no eixo y
    Li= sqrt((Lxi^2)+(Lyi^2));  % Comprimento do elemento
    Lvec(i)= Li;                % Vetor de comprimentos
    
    C= Lxi/Li;                  % Angulo cosseno do elemento
    S= Lyi/Li;                  % Angulo seno do elemento
    Cvec(i)= C;                 % Vetor dos cossenos C
    Svec(i)= S;                 % Vetor dos senos C
    
    Ke= [C^2 C*S -C^2 -C*S;     % Matriz de rigidez do elemento
         C*S S^2 -C*S -S^2;
        -C^2 -C*S C^2 C*S;
        -C*S -S^2 C*S S^2]*eal; 
    
%     Me= [2*C^2 2*C*S -C^2 -C*S;	% Matriz de massa consistente
%          2*C*S 2*S^2 -C*S -S^2;
%         -C^2 -C*S 2*C^2 2*C*S;
%         -C*S -S^2 2*C*S 2*S^2]*(rho*A*Li)/6; 
%     
%     Ml= [C^2 C*S  0   0;        % Matriz de massa concentrada
%          C*S S^2  0   0;
%           0   0  C^2 C*S;
%           0   0  C*S S^2]*(rho*A*Li)/2; 
    
    u1= gdls(no1,1);            % Restricao em x do elemento 1
    v1= gdls(no1,2);            % Restricao em x do elemento 1
    u2= gdls(no2,1);            % Restricao em x do elemento 1
    v2= gdls(no2,2);            % Restricao em x do elemento 1
    
    gdlEli= [u1 v1 u2 v2];      % GDL do elemento i
    
    % Montando matriz de rigidez global para os elementos livres
    for j=1:length(gdlEli)
        for k=1:length(gdlEli)
            if (gdlEli(j)~=0 && gdlEli(k)~=0)
                Keg(gdlEli(j),gdlEli(k))= Keg(gdlEli(j),gdlEli(k)) + Ke(j,k);
            end
        end
    end
end

Fvec= zeros(Ngdl,1);
for i=1:size(F,1)
    noF= F(i,2);                % No da forca aplicada
    gdlF= F(i,3);               % GDL da forca aplicada
    gdlFG= gdls(noF,gdlF);
    Fvec(gdlFG,1)= F(i,4);
end

%% Resolvendo o sistema

desloc= Keg\Fvec;               % Resolvendoos deslocamentos

nosdesloc= nos;                 % Iniciando vetor de nos deslocados
for i=1:size(gdls,1)
    for j=1:size(gdls,2)
        if gdls(i,j)~=0
            nosdesloc(i,j+1)= nosdesloc(i,j+1) + desloc(gdls(i,j),1);
        end
    end
end

Lvec2= zeros(1,size(el,1));     % Comprimento apos deslocamento
for i=1:size(el,1)
    no1= el(i,2);           
    no2= el(i,3);
    x1= nosdesloc(no1,2);       % Posicao em x do no 1 do elemento
    x2= nosdesloc(no2,2);       % Posicao em x do no 2 do elemento
    y1= nosdesloc(no1,3);       % Posicao em y do no 1 do elemento
    y2= nosdesloc(no2,3);       % Posicao em y do no 2 do elemento
    
%     Lxi= abs(x2-x1);                 % Comprimento do elemento no eixo x
%     Lyi= abs(y2-y1);                 % Comprimento do elemento no eixo y
    
    Lxi= x2-x1;                 % Comprimento do elemento no eixo x
    Lyi= y2-y1;                 % Comprimento do elemento no eixo y
    Li= sqrt((Lxi^2)+(Lyi^2));  % Comprimento do elemento
    Lvec2(i)= Li;               % Vetor de comprimentos
    
    C= Lxi/Li;                  % Angulo cosseno do elemento
    S= Lyi/Li;                  % Angulo seno do elemento
    Cvec2(i)= C;                % Vetor dos cossenos C
    Svec2(i)= S;                % Vetor dos senos C
end

dL= Lvec2-Lvec;                 % Diferenca de comprimento, ou deslocamento
Kgf= zeros(Ngdl,Ngdl);          % Matriz de rigidez simetrica em coord Global

% A= 2.2361;
% E= 2.2361;

for i=1:length(dL)
%     Fn(i)= A*E*dL(i)/Lvec(i);  
%     C= Cvec2(i);
%     S= Svec2(i);

    Fn(i)= dL(i);               % A*E*dL(i)/Lvec(i) - VERIFICAR!!!
    C= Cvec(i);
    S= Svec(i);
  
    Kg=[S^2 -C*S -S^2  C*S;     % Matriz de rigidez em coord global
       -C*S  C^2  C*S -C^2;
       -S^2  C*S  S^2 -C*S;
        C*S -C^2 -C*S  C^2];
    
    Kg= Fn(i)/Lvec(i)*Kg;
    
    no1= el(i,2);
    no2= el(i,3);
    u1= gdls(no1,1);            % Restricao em x do elemento 1
    v1= gdls(no1,2);            % Restricao em x do elemento 1
    u2= gdls(no2,1);            % Restricao em x do elemento 1
    v2= gdls(no2,2);            % Restricao em x do elemento 1
    gdlEli2= [u1 v1 u2 v2];     % GDL do elemento i
    
    % Montando matriz de rigidez em coordenadas globais
    for j=1:length(gdlEli2)
        for k=1:length(gdlEli2)
            if (gdlEli2(j)~=0 && gdlEli2(k)~=0)
                Kgf(gdlEli2(j),gdlEli2(k))= Kgf(gdlEli2(j),gdlEli2(k)) + Kg(j,k);
            end
        end
    end    
end
%% Analise modos e resultados

[fmodal,freqmodo]= eig(Keg,Kgf);% Modos

fmodal1= fmodal(:,1);

nosDeslf1= nos;                  % Iniciando vetor de nos deslocados
for i=1:size(gdls,1)
    for j=1:size(gdls,2)
        if gdls(i,j)~=0
            nosDeslf1(i,j+1)= nosDeslf1(i,j+1) + fmodal1(gdls(i,j),1)/2;
        end
    end
end
% Plot da trelica no modo
for i=1:size(el,1)
    x1=nosDeslf1(el(i,2),2);
    x2=nosDeslf1(el(i,3),2);
    y1=nosDeslf1(el(i,2),3);
    y2=nosDeslf1(el(i,3),3);
    figure(1)
    line([x1,x2],[y1,y2], 'Color', 'b', 'linestyle', '-', 'linewidth', 1.5); 
end

fmodal2= fmodal(:,2);
nosDeslf2= nos;                  % Iniciando vetor de nos deslocados
for i=1:size(gdls,1)
    for j=1:size(gdls,2)
        if gdls(i,j)~=0
            nosDeslf2(i,j+1)= nosDeslf2(i,j+1) + fmodal2(gdls(i,j),1)/2;
        end
    end
end

for i=1:size(el,1)
    x1=nosDeslf2(el(i,2),2);
    x2=nosDeslf2(el(i,3),2);
    y1=nosDeslf2(el(i,2),3);
    y2=nosDeslf2(el(i,3),3);
    figure(1)
    line([x1,x2],[y1,y2], 'Color', 'r', 'linestyle', '-', 'linewidth', 1.5);
end
axis equal
title('Truss buckling modes')
xlabel('x')
ylabel('y')
legend('Undeformed','Mode 1','Mode 2', 'location', 'northwest')

%%

for i=1:size(el,1)
    x1=nos(el(i,2),2);
    x2=nos(el(i,3),2);
    y1=nos(el(i,2),3);
    y2=nos(el(i,3),3);
    figure(2)
    line([x1,x2],[y1,y2], 'Color', 'k', 'linestyle', '--', 'linewidth', 1.5); 
    hold on; grid on
end
for i=1:size(el,1)
    x1=nosdesloc(el(i,2),2);
    x2=nosdesloc(el(i,3),2);
    y1=nosdesloc(el(i,2),3);
    y2=nosdesloc(el(i,3),3);
    figure(2)
    line([x1,x2],[y1,y2], 'Color', 'r', 'linestyle', '-', 'linewidth', 1.5);
end
axis equal
title('Truss deformation')
xlabel('x')
ylabel('y')

