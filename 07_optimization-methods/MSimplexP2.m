% Método Simplex
%%
clear all; close all; clc;

%% Problema 2 Lista 8
A = [   1	1	1	0	0	0	0	0	0	1	0	0 ;...
        0	0	0	1	1	1	0	0	0	0	1	0 ;...
        0	0	0	0	0	0	1	1	1	0	0	1 ;...
        1	0	0	1	0	0	1	0	0	0	0	0 ;...
        0	1	0	0	1	0	0	1	0	0	0	0 ;...
        0	0	1	0	0	1	0	0	1	0	0	0 ];
    
b = [ 2000 ; 3000 ; 1500 ; 2000 ; 2000 ; 1000 ];
c = -[25 20 30 30 25 25 20 15 23 0 0 0];
cini = -[25 20 30 30 25 25 20 15 23 0 0 0];
I = [7 8 9 10 11 12];    
T = 1:1:size(A,2);
J = setdiff(T,I);
    
% Z= -x1+2*x2;
mA= length(A(:,1));                 %Dimensão m de A
nA= length(A(1,:));                 %Dimensão n de A
% x= [x1, x2, x3, x4, x5];          %Vetor de variáveis (Apenas representativo)


nI= length(I);                      %Dimensão n de I
nJ= length(J);                      %Dimensão n de J

CheckDim= nI+nJ;                    %Check do tamanho das dimensões I e J
if CheckDim~=nA
    errmsg=(sprintf('Erro, dim(I)+dim(J) ~= n. \nValor resultante = %d  n= %d.', CheckDim, nA));
    error(errmsg);
end

%Início primeira parte algoritmo
%__________________________________________________________________________

Ai= A(:,I);                         %Separando submatriz Ai de A
Aj= A(:,J);                         %Separando submatriz Aj de A

nAi= length(A(:,1));                %Dimensão n de Ai

ID= eye(nAi);
% ID= [1 0 0;...                      %Matriz Identidade
%      0 1 0;...
%      0 0 1];

for i= 1:length(ID)
    for j= 1:length(ID)
        if Ai(i,j)==ID(i,j)
        	compID(i,j)= 1;         %Comparação Ai x ID (igual = 1, deferente = 0)
        else                        
            compID(i,j)= 0;
        end
    end
end

cont=0;                             %Iniciando contador
cond=0;                             %Iniciando teste de condição

for k= 1:nAi
    indZero= find(compID(:,k)==0);          %Procurando indices com zero
    
    if isempty(indZero)==0                  %Testa indice nulo
        cond= 1;                            %Alterando condicional
        cont= cont+1;                       %Incremento do contador
        J= [J I(k)];                        %Alterando vet de var não básicas
        I(k)= nI+I(k);                      %Alterando vet de var básicas
        Aj= [Aj Ai(:,k)];                   %Alterando Aj
        Ai(:,k)= ID(:,k);                   %Igualando Ai com ID
    end
end

contpt1= 0;
ctemp= c;                           %Variável c temporário
if cond==1                          %Minimizando Phi se cond = 1
    c2= zeros(1,CheckDim);          %Zerando c
    ctemp= [ctemp zeros(1,cont)];
    c2(I)= -1;                      %Adicionando -1 em c nas variáveis I
    A= [Aj Ai];                     %Reescrevendo a matriz A
    cond2= 1;                       
    I= sort(I);                     %Ordenando o vetor I
    base0= zeros(1,length(c2));      
	
    for i=1:nI
        ind1= find(A(:,I(i))==1);
        piv= A(ind1,I(i))/c2(I(i));
        c2= c2-piv*A(ind1,:);
    end
    
    while cond2==1
        Ai= A(:,I);                 %Separando A de I
        Aj= A(:,J);                 %Separando A de J
        ci= c2(I);                  %Separando c de I
        cj= c2(J);                  %Separando c de J
        
        varreaj= ci/Ai;             %Variavel de reajuste
        cj= cj-varreaj*Aj;          %Reajuste de cj
        
        pcj= cj(cj>0);              %Separando variáveis positivas de cj
        
        if isempty(pcj)==0
            [~,indmaxcj]= max(cj);          %Indice do máximo de cj
            indJ= J(indmaxcj);              %Indice de J relacionado a indmaxcj
            base0(indJ)= 1;
            valAindJ= A(:,indJ);            %Fator de correção 
            alfa= b./valAindJ;              %Aplicando fator em b
            [alfac,alfaind]= sort(alfa);    %Organizando em ordem crescente
            alfaind= alfaind(alfac>=0);     
            
            for ii=1:length(alfaind)
                valalfa= alfaind(ii);
                if base0(I(valalfa))==0
                    valIalfa= I(valalfa);
                    break;
                end                
            end
            
            IAb= [I' A b];         %Matriz completa I A b
            
            %Pivoteamento--------------------------------------------------
            
            indvalI= find(IAb(:,1)==valIalfa);              %Linha para substituir
            Ab= IAb(:,2:end);                               %Unindo A e b
            Ab(indvalI,:)= Ab(indvalI,:)/Ab(indvalI,indJ);  %Pivoteando linha
            
            nAb= size(Ab,1);
            for jj= 1:nAb           %Operando as outras linhas
                if Ab(jj,indJ)~=0 && jj~=indvalI
                    Ab(jj,:)= Ab(jj,:)-Ab(jj,indJ)*Ab(indvalI,:);
                end
            end
            c2= c2-c2(indJ)*Ab(indvalI,1:end-1);   %
            A= Ab(:,1:end-1);       %Separando A de Ab
            b= Ab(:,end);           %Separando b de Ab
            
            ctemp= ctemp-ctemp(indJ)*Ab(valalfa,1:end-1);
            
            %--------------------------------------------------------------
%             Ap1= A
%             bp1= b
            
            indtemp= I==valIalfa;   %Trocando as bases
            I(indtemp)= indJ;
            indtemp2= J==indJ;
            J(indtemp2)= valIalfa;
        else
            cond2=0;
            Ai= A(:,I);
            varreaj= ci/Ai;
            Z= varreaj*b;
            result= inv(Ai)*b;
        end
        contpt1= contpt1+1;
    end
%     App1= A
%     bpp1= b
%     b
    
    A= A(:,c2==0);                  %Separando variáveis básicas
    
%     A2= A
%     b
    
    ctemp= ctemp(:,c2==0);
    
    aux1=0;
    aux2=1;
    N= length(A(:,1));
    
    for l=1:N
        cond3=1;
        while aux2<=length(A(1,:)) && cond3==1
            tlog= A(:,aux2)==ID(:,l);           %Comparando A com mat Identidade
            tlog= tlog*1;                       %Convertendo para número
            if mean(tlog)==1
                aux1= aux1+1;
                I(aux1)= aux2;
                cond3=0;
            end
            aux2= aux2+1;
        end
        aux2= 1;
    end
    
    T= 1:1:size(A,2);
    J= setdiff(T,I);
end

% A3= A

%Fim primeira parte algoritmo
%__________________________________________________________________________

c=ctemp;

cpts= 1;
base1= zeros(1,length(c));
contAB= 0;
contpt2= 0;

while cpts==1
	Ai= A(:,I);                 %Separando A de I
    Aj= A(:,J);                 %Separando A de J
    ci= c(I);                   %Separando c de I
	cj= c(J);                   %Separando c de J
    
%     I
%     J
%     PreAp2= A
% 	  Prebp2= b
        
	varreaj= ci/Ai;             %Variavel de reajuste
    cj= cj-varreaj*Aj;          %Reajuste de cj

    pcj2= c(c>0);               %Separando variáveis positivas de cj
        
	if isempty(pcj2)==0
        [~,indmaxcj]= max(cj);              %Indice do máximo de cj
        indJ= J(indmaxcj);                  %Indice de J relacionado a indmaxcj
        base1(indJ)= 1;
        valAindJ= A(:,indJ);                %Fator de correção 
        alfa2= b./valAindJ;                 %Aplicando fator em b
        [alfac2,alfaind2]= sort(alfa2);     %Organizando em ordem crescente
        alfaind2= alfaind2(alfac2>=0);     
            
        for ii=1:length(alfaind2)
            valalfa2= alfaind2(ii);
            if base1(I(valalfa2))==0
            	valIalfa2= I(valalfa2);
                break;
            end                
        end
            
        IAb2= [I' A b];         %Matriz completa I A b
        
        %Pivoteamento2-----------------------------------------------------
            
        indvalI= find(IAb2(:,1)==valIalfa2);    %Linha para substituir
        
%         if isempty(indvalI)==1
%             cpts=0;
%             break
%         end
        
        Ab= IAb2(:,2:end);                      %Unindo A e b
        
%         if Ab(indvalI,indJ)~=0
            Ab(indvalI,:)= Ab(indvalI,:)./Ab(indvalI,indJ); %Pivoteando linha
%         else
%             Ab(indvalI,:)= Ab(indvalI,:)./1; %Pivoteando linha
%         end
        
        nAb= size(Ab,1);
        for jj= 1:nAb                           %Operando as outras linhas
            if Ab(jj,indJ)~=0 && jj~=indvalI
                Ab(jj,:)= Ab(jj,:)-Ab(jj,indJ)*Ab(indvalI,:);
            end
        end

%         c= c(1:4)-c(indJ)*Ab(indvalI,1:end-1);   %        
        c= c-c(indJ)*Ab(indvalI,1:end-1);   %        

        A= Ab(:,1:end-1);           %Separando A de Ab
        b= Ab(:,end);               %Separando b de Ab
        
%         ctemp= ctemp-ctemp(indJ)*Ab(valalfa,1:end-1);
        ctemp= ctemp-ctemp(indJ)*Ab(valalfa2,1:end-1);
        
        %------------------------------------------------------------------
%         Ap2= A
%         bp2= b
        
        indtemp2= I==valIalfa2;     %Trocando as bases
        I(indtemp2)= indJ;
        indtemp21= J==indJ;
        J(indtemp21)= valIalfa2;
        
        A
%         c
%         b
%         valIalfa2
    else
        cpts=0;
        Ai= A(:,I);
        
        ci= cini(I);
        
        varreaj2= ci/Ai;
        Z= varreaj2*b;
        
        result= inv(Ai)*b;
    end
    contpt2= contpt2+1;
    
%     A
%     b
%     ctemp
    
end
% App2= A
% bpp2= b

resultFinal= zeros(1,length(A(1,:)));
resultFinal(I)= result;

fprintf('Função: z= %.4f \n',Z);
for o=1:length(resultFinal)
    fprintf('x%d= %.4f \n',o,resultFinal(o));
end

