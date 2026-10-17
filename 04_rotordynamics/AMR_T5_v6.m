% Analise de maquinas rotativas - Exercicio 5 - Modelo de eixo por elementos finitos e mancais com coeficientes dinamicos
% Nicolas Londero Gass - nicolas.gass@eac.ufsm.br - 2022/2°S
clear all; close all; clc
warning off;

%% Dados iniciais

% Dados do eixo
De= 13*10^(-3);                                     % Diametro do eixo [m]
Ce= 747*10^(-3);                                    % Comprimento do eixo [m]
Ae= pi*(De/2)^2;                                    % Area do transversal do eixo [m²]
Ie= (pi*(De/2)^4)/4;                                % Momento de inercia eixo [m^4]

% Dados do mancal
Dm= 30*10^-3;                                       % Diametro do mancal [m]
R= Dm/2;                                            % Raio do mancal [m]
Cm= 20*10^-3;                                       % Comprimento mancal [m]
Am= pi*(R)^2;                                       % Area do transversal do mancal [m²]
frd= 90*10^(-6);                                    % Folga Radial [m]
eta= 0.051;                                         % Viscosidade absoluta (30°C)
Im= (pi*(R)^4)/4;                                   % Momento de inercia mancal [m^4]

% Dados do disco
Dd= 90*10^(-3);                                     % Diametro do disco [m]
Cd= 47*10^(-3);                                     % Comprimento do disco [m]
Ad= pi*(Dd/2)^2;                                    % Area do disco [m²]
Vd= Ad*Cd;                                          % Volume do disco [m^3]
rho= 7850;                                          % Densidade do aco [kg/m³]
md= Vd*rho;                                         % Massa do disco [kg];
Ip= (md*(Dd^2))/8;                                  % Momento de inercia polar do disco [kg.m²]
Id= Ip/2 + (md*Cd^2)/12;                         % Momento de inercia diametral do disco [kg.m²]
a= 145*10^(-3);                                     % Distancia do disco para o mancal 1 no eixo
b= Ce-a;                                            % Distancia do disco para o mancal 2 no eixo

% Dados gerais
ME= 200*10^9;                                       % Modulo de elasticidade do aco [Pa]
ps= 0.3;                                            % Coef Poisson [-]
ke= (48*ME*Ie)/(Ce^3);                              % Constante elastica do eixo [N/m]
c= 10^(-4);                                         % Coef de amortecimento [N.s/m]
me= 3.7*10^(-5);                                    % Desbalanceamento [kg.m]
g= 9.81;                                            % Gravidade [m/s²]

% Analise estatica 
Fp= md*g;                                           % Forca peso aplicado pela massa do disco
Fm1= Fp*(b/Ce);                                     % Distribuicao da forca peso no mancal 1
Fm2= Fp*(a/Ce);                                     % Distribuicao da forca peso no mancal 2
F0= [Fm1, Fm2];                                     % Forcas aplicadas aos mancais    

omgrd= 10:5:1500;                                   % Vetor omega em radianos
omghz= omgrd/(2*pi);                                % Vetor omega em Hertz

valest= [0.6 0.1; 0.6 0.1];                         % Valores estimados de epsilon e phi, tirado do slide 39

n= size(omgrd,2);                                   

% Inicializando vetores
Phi= zeros(1,n); phi= zeros(1,n); Sm22= zeros(1,n); 
Sm2= zeros(1,n); g22= zeros(1,n); g21= zeros(1,n); 
g12= zeros(1,n); g11= zeros(1,n); b22= zeros(1,n); 
b21= zeros(1,n); b12= zeros(1,n); b11= zeros(1,n); 
k11= zeros(1,n); k12= zeros(1,n); k21= zeros(1,n); 
k22= zeros(1,n); d11= zeros(1,n); d12= zeros(1,n); 
d21= zeros(1,n); d22= zeros(1,n); E3= zeros(1,n); 
Km1= zeros(4,4,n); Cm1= zeros(4,4,n);
Km2= zeros(4,4,n); Cm2= zeros(4,4,n);
ey= zeros(2,n); ez= zeros(2,n);

% Waitbar ------------------------------------------
steps = length(omgrd)*length(F0);
f = waitbar(0.0001,sprintf('Passo %d de %d - %5.2f %%',0,steps,(0.0001)*100));
setappdata(f,'canceling',0);
cntwb= 1;                                                                       % Contador do waitbar
% --------------------------------------------------

tic
for i=1:length(F0)
    valest2= valest(i,:)';
    for j=1:length(omgrd)
    
        if getappdata(f,'canceling')                                            % Cancelamento waitbar
            break
        end

        K0= (eta*omgrd(j)*R*(Cm^3))/(frd^2);        % Fator de forca do mancal curto
        [E3(i,j),Phi(i,j)]= MNR_bearing(K0, F0(i), valest2);      % Epsilon e Phi por Newton-Raphson (Jacobiano analitico)

        valest2= [E3(i,j); Phi(i,j)];                                           % Variaveis finais
        phi(i,j)= Phi(i,j)*(180/pi);                                            % phi em graus

        Fn2= (eta*Cm^3*omgrd(j)*(R))/(2*frd^2);                                 % Forca de fricao para omghz(i)

        Sm2(i,j)= F0(i)/Fn2;                                                    % Sommerfeld modificado
        Sm22(i,j)= (pi/2)*(E3(i,j)/(1-E3(i,j)^2)^2)*sqrt(1-E3(i,j)^2+((4/pi)*E3(i,j))^2); % Sommerfeld modificado

        % Calculos de Kik e Cik e variaveis auxiliares
        A= 4/(pi^2+(16-pi^2)*E3(i,j)^2)^(3/2);                                      

        g11(i,j)= (2*pi^2+(16-pi^2)*E3(i,j)^2)*A;                                                                   % Gamma11
        g12(i,j)= (pi/4)* ((pi^2-(2*pi^2)*E3(i,j)^2-(16-pi^2)*E3(i,j)^4)/ (E3(i,j)*sqrt(1-E3(i,j)^2)))*A;           % Gamma12
        g21(i,j)= (-pi/4)* ((pi^2+(32+pi^2)*(E3(i,j)^2)+(32-2*pi^2)*(E3(i,j)^4))/ (E3(i,j)*sqrt(1-E3(i,j)^2)))*A;   % Gamma21
        g22(i,j)= (((pi^2)+(32+(pi^2))*(E3(i,j)^2)+ (32-2*(pi^2))*(E3(i,j)^4))/ (1-(E3(i,j)^2)))*A;                  % Gamma22 (denominador 1-eps^2)
                
        b11(i,j)= (pi/2)*((sqrt(1-E3(i,j)^2))/E3(i,j))* (pi^2+(2*pi^2-16)*E3(i,j)^2)*A;                             % Beta11
        b12(i,j)= -(2*pi^2+(4*pi^2-32)*E3(i,j)^2)*A;                                                                % Beta12
        b21(i,j)= b12(i,j);                                                                                         % Beta21
        b22(i,j)= (pi/2)*((pi^2+(48-2*pi^2)*E3(i,j)^2+(pi^2)*E3(i,j)^4)/ (E3(i,j)*sqrt(1-E3(i,j)^2)))*A;            % Beta22

        k11(i,j)= g11(i,j)*(F0(i)/frd);                     % k11
        k12(i,j)= g12(i,j)*(F0(i)/frd);                     % k12
        k21(i,j)= g21(i,j)*(F0(i)/frd);                     % k21
        k22(i,j)= g22(i,j)*(F0(i)/frd);                     % k22

        d11(i,j)= b11(i,j)*(F0(i)/(frd*omgrd(j)));          % d11
        d12(i,j)= b12(i,j)*(F0(i)/(frd*omgrd(j)));          % d12
        d21(i,j)= b21(i,j)*(F0(i)/(frd*omgrd(j)));          % d21
        d22(i,j)= b22(i,j)*(F0(i)/(frd*omgrd(j)));          % d22
        
        % Matriz K e C de cada mancal
        if i==1
            Km1(:,:,j)=[k11(i,j) 0 k12(i,j) 0;
                           0     0    0     0;
                        k21(i,j) 0 k22(i,j) 0;
                           0     0    0     0];
                   
            Cm1(:,:,j)=[d11(i,j) 0 d12(i,j) 0;
                           0     0    0     0;
                        d21(i,j) 0 d22(i,j) 0;
                           0     0    0     0];
                   
        else
            Km2(:,:,j)=[k11(i,j) 0 k12(i,j) 0;
                           0     0    0     0;
                        k21(i,j) 0 k22(i,j) 0;
                           0     0    0     0];
                   
            Cm2(:,:,j)=[d11(i,j) 0 d12(i,j) 0;
                           0     0    0     0;
                        d21(i,j) 0 d22(i,j) 0;
                           0     0    0     0];
            
        end
        
        % Update waitbar e mensagem
        waitbar(cntwb/steps,f,sprintf('Passo %d de %d - %4.2f %%',cntwb,steps,(cntwb/steps)*100));
        cntwb=cntwb+1;
        %--------------------------
    end
end
tt(1)= toc;
delete(f)

for k=1:2
    ey(k,:)= (E3(k,:)*frd/2).*sind(phi(k,:));
    ez(k,:)= -(E3(k,:)*frd/2).*cosd(phi(k,:));    
end
for k=1:2
    ey2(k,:)= (E3(k,:)).*sind(phi(k,:));
    ez2(k,:)= -(E3(k,:)).*cosd(phi(k,:));    
end

%% Analise em elementos finitos

caso=3;
%%% Definir discretizacao desejada %%%
% Caso 1: 15 nos
% Caso 2: 18 nos
% Caso 3: 30 nos

% n elementos -> x
% no de posição do mancal e do disco
if caso==1
    % Coordenadas fisicas dos nos (Configuracao 1)
    %%%       ----|  M  |  M  |----.----|   D  |   D  |-------.-------.-------.-------.-------.-------|   M  |   M  |------
    %            -10    0    10  72.5  121.5  145   168.5   217.5    290    362.5    435    507.5    570    580    590                                
    %                   x                                                                                    x    
    %         0                                                                                                                 747
    %%%       -----|  M  |  M  |-------.-------|   D  |   D  |-------.-------.-------.-------.-------.-------|   M  |   M  |-----
    %         0   73.5  83.5  93.5    156     195     x     242     301    373.5    446    518.5    591    653.5  663.5  673.5  747
    %   >>>   0   73.5  83.5  93.5    156      x    228.5    x      301    373.5    446    518.5    591    653.5  663.5  673.5  747
    cordf= [0 73.5 83.5 93.5 156 228.5 301 373.5 446 518.5 591 653.5 663.5 673.5 747]*0.001;
        
    % Dimensoes dos elementos (Configuracao 1)
    Ael= [Ae, Am, Am, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Am, Am, Ae];
    Del= [De, Dm, Dm, De, De, De, De, De, De, De, De, De, De, Dm, Dm, De];
    Iel= [Ie, Im, Im, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Im, Im, Ie];
    
    % Definição da posição dos mancais e do disco (Configuracao 1)
    LMc1= 3;            % Linha da matriz onde se encontra o centro do mancal 1 (elemento da matriz de incidencia)
    LMc2= 13;           % Linha da matriz onde se encontra o centro do mancal 2 (elemento da matriz de incidencia)
    LDc= 6;             % Linha da matriz onde se encontra o centro do disco (elemento da matriz de incidencia)
        
elseif caso==2
    % Coordenadas fisicas dos nos (Configuracao 2)
    %%%     ----|  M  |  M  |----------|   D  |   D  |-----------------------------------------------|   M  |   M  |------
    cordf= [0 51.25 61.25 71.25 127 182.75 206.25 229.75 285.5 341.25 397 452.75 508.5 564.25 620 675.75 685.75 695.75 747]*0.001;
    
    % Dimensoes dos elementos (Configuracao 2)
    Ael= [Ae, Am, Am, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Am, Am, Ae];
    Del= [De, Dm, Dm, De, De, De, De, De, De, De, De, De, De, De, De, Dm, Dm, De];
    Iel= [Ie, Im, Im, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Im, Im, Ie];
    
    % Definição da posição dos mancais e do disco (Configuracao 2)
    LMc1= 3;            % Linha da matriz onde se encontra o centro do mancal 1 (elemento da matriz de incidencia)
    LMc2= 17;           % Linha da matriz onde se encontra o centro do mancal 2 (elemento da matriz de incidencia)
    LDc= 7;             % Linha da matriz onde se encontra o centro do disco (elemento da matriz de incidencia)
    
elseif caso==3
    % Coordenadas fisicas dos nos (Configuracao 32)
    %%%     ----.----|  M  |  M  |------.------.------.------|   D  |   D  |-------.-------.-------.-------.-------.-------.
    %       0 36.75 73.5  83.5  93.5  119.75  156   192.25   x    228.5    x    264.75    301    337.25  373.5  409.75    446   
    %            482.25  518.5   554.75   591    627.25  653.5  663.5  673.5  710.25 747                                             
    %       -------.-------.-------.-------.-------.-------|   M  |   M  |------.------              
    cordf= [0 36.75 73.5 83.5 93.5 119.75 156 192.25 228.5 264.75 301 337.25 373.5 409.75 446 482.25 518.5 554.75 591 627.25...
            653.5 663.5 673.5 710.25 747]*0.001;
    
    % Dimensoes dos elementos (Configuracao 3)
    Ael= [Ae, Ae, Am, Am, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Ae, Am, Am, Ae, Ae];
    Del= [De, De, Dm, Dm, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, De, Dm, Dm, De, De];
    Iel= [Ie, Ie, Im, Im, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Ie, Im, Im, Ie, Ie];
    
    % Definição da posição dos mancais e do disco (Configuracao 3)
    LMc1= 4;            % Linha da matriz onde se encontra o centro do mancal 1 (elemento da matriz de incidencia)
    LMc2= 22;           % Linha da matriz onde se encontra o centro do mancal 2 (elemento da matriz de incidencia)
    LDc= 9;            % Linha da matriz onde se encontra o centro do disco (elemento da matriz de incidencia)
        
end

% Elementos
ne= length(cordf)-1;                                % Numero de elementos
nn= ne+1;                                           % Numero de nos
ngdl= 4;                                            % Numero de graus de liberdade por elemento
ngdlt= nn*4;                                        % Numero de graus de liberdade total
ggdl= [1:nn, nn+1:ngdlt];                           % Numero de graus de liberdade globais
vecel= 1:ngdl;                                      % Vetor de gdl dos elementos na mat de incidencia
% ^^ Vetor de inserção dos elementos dos mancais e disco (matriz 4x4 => 4 el)

% Inicializando vetores das matrizes globais, matriz de incidencia e deslocamentos
Fg= zeros(ngdlt,1); incel= zeros(ne, 8);
Mg= zeros(ngdlt); Cg= zeros(ngdlt);
Kg= zeros(ngdlt); Gg= zeros(ngdlt); 
X= zeros(ngdlt, length(omgrd));

% Waitbar ------------------------------------------
steps = length(1:ne)+length(omgrd);
f = waitbar(0.0001,sprintf('Passo %d de %d - %5.2f %%',0,steps,(0.0001)*100));
setappdata(f,'canceling',0);
cntwb= 1;                                           % Contador do waitbar
% --------------------------------------------------
    
% Montando matrizes globais de rigidez, massa, amortecimento e giroscopica
for el=1:ne
    
    Le= cordf(el+1)-cordf(el);                      % Comprimento de cada elemento (entre nos)
    
    % Matriz de rigidez do elemento de viga com 8GDL
    Kel= [  12     6*Le     0        0       -12    6*Le     0        0;
          6*Le    4*Le^2    0        0     -6*Le   2*Le^2    0        0;
            0        0      12     -6*Le      0       0     -12    -6*Le;
            0        0    -6*Le    4*Le^2     0       0     6*Le   2*Le^2;
           -12    -6*Le     0        0        12   -6*Le     0        0;
          6*Le    2*Le^2    0        0      -6*Le  4*Le^2    0        0;
            0        0     -12      6*Le      0       0      12      6*Le;
            0        0    -6*Le    2*Le^2     0       0     6*Le    4*Le^2] *((ME*Iel(el))/Le^3);
    
    % Matriz de massa do elemento de viga com 8GDL
	Mel= [ 156      22*Le     0         0        54    -13*Le     0        0;
          22*Le    4*Le^2     0         0      13*Le  -3*Le^2     0        0;
            0         0      156     -22*Le      0        0       54     13*Le;
            0         0    -22*Le     4*Le^2     0        0     -13*Le  -3*Le^2;
            54      13*Le     0         0       156    -22*Le     0        0;
         -13*Le   -3*Le^2     0         0     -22*Le  4*Le^2      0        0;
            0         0       54     -13*Le      0        0      156     22*Le;
            0         0     13*Le    -3*Le^2     0        0      22*Le   4*Le^2] *((rho*Ael(el)*Le)/420);
                
	% Matriz giroscopica de viga com 8GDL
    Gel= [ 0      0      36   -3*Le    0      0      -36   -3*Le;
           0      0    3*Le  -4*Le^2   0      0     -3*Le   Le^2;
          -36   -3*Le    0      0      36   -3*Le     0       0;
          3*Le  4*Le^2   0      0    -3*Le  -Le^2     0       0;
           0      0     -36    3*Le    0      0       36     3*Le;
           0      0     3*Le   Le^2    0      0     -3*Le  -4*Le^2;
           36    3*Le    0      0     -36    3*Le     0       0
          3*Le  -Le^2    0      0    -3*Le  4*Le^2    0       0] *(rho*Ael(el)*Del(el)^2)/(240*Le);
    
    % Amortecimento do elemento
    alp= 0;
    beta= 1*10^(-4);
    Cel= alp*Mel+beta*Kel;
    
    % Vetor de incidencia do elemento (Graus de liberdade de cada elemento)
    incel(el,:)= [ggdl(4*el-3:4*el), ggdl(4*(el+1)-3:4*(el+1))];             % Estudar !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
	% Montando o elemento na matriz global
    Kg(incel(el,:),incel(el,:))= Kg(incel(el,:),incel(el,:)) + Kel;
    Mg(incel(el,:),incel(el,:))= Mg(incel(el,:),incel(el,:)) + Mel;
    Gg(incel(el,:),incel(el,:))= Gg(incel(el,:),incel(el,:)) + Gel;
    Cg(incel(el,:),incel(el,:))= Cg(incel(el,:),incel(el,:)) + Cel;	
    
    % Update waitbar e mensagem
	waitbar(cntwb/steps,f,sprintf('Passo %d de %d - %4.2f %%',cntwb,steps,(cntwb/steps)*100));
    cntwb= cntwb+1;
    %--------------------------
end

% Adicionando os efeitos nas matrizes globais
Md= [md 0  0  0                                     % Efeito de inercia e massa do disco
	 0  Id 0  0                                     
	 0  0  md 0
	 0  0  0  Id];
Mg(incel(LDc,vecel),incel(LDc,vecel))= Mg(incel(LDc,vecel),incel(LDc,vecel)) + Md;
    
Gd= [0  0  0  0                                     % Efeito giroscopico do disco
	 0  0  0  Ip
	 0  0  0  0
	 0 -Ip 0  0];
Gg(incel(LDc,vecel),incel(LDc,vecel))= Gg(incel(LDc,vecel),incel(LDc,vecel)) + Gd;

% Incluindo fatores adicionais nas matrizes globais
for k=1:length(omgrd)
    Kg2= Kg; Mg2= Mg;
    Gg2= Gg; Cg2= Cg;
    Fg= zeros(ngdlt,1);
    s= 1i*omgrd(k);                                 % Variavel de Laplace (TL)
    
    % Vetor de forcas aplicadas (desbalanceamento girante: y e z defasados de 90 graus)
    % O peso (Fp) e uma carga estatica, ja considerada no equilibrio dos mancais (F0),
    % e nao entra na resposta harmonica.
    Fg(incel(LDc,1),1)= me*omgrd(k)^2;
    Fg(incel(LDc,3),1)= -1i*me*omgrd(k)^2;
    
    % Adicionando rigidez do mancal 1
    Kg2(incel(LMc1, vecel), incel(LMc1, vecel))= Kg2(incel(LMc1, vecel), incel(LMc1, vecel)) + Km1(:,:,k);
	
    % Adicionando rigidez do mancal 2
    Kg2(incel(LMc2, vecel), incel(LMc2, vecel))= Kg2(incel(LMc2, vecel), incel(LMc2, vecel)) + Km2(:,:,k);
    
    % Adicionando amortecimento do mancal 1
    Cg2(incel(LMc1, vecel), incel(LMc1, vecel))= Cg2(incel(LMc1, vecel), incel(LMc1, vecel)) + Cm1(:,:,k);
    
    % Adicionando amortecimento do mancal 2
    Cg2(incel(LMc2, vecel), incel(LMc2, vecel))= Cg2(incel(LMc2, vecel), incel(LMc2, vecel)) + Cm2(:,:,k);
	
    % Equacao do deslocamento
    X(:,k)= ((Mg2*s^2)+(Cg2+(omgrd(k)*Gg2))*s+Kg2)\Fg;
    %%% X -> Coluna=> Frequência; Linha=> GDL de cada elemento
    %%% 4gdl por elemento - primeiro GDL - deslocamento em y
    %                     - segundo GDL - rotacao em y
    %                     - terceiro GDL - deslocamento em z
    %                     - quarto GDL - rotação em z
    %
    %%% Deslocamento em y - primeiro elemento (e1y) linha 1, e2y- l 5,
    %%% e3y- l 9, e4y- l 13, e5y- l 17, e6y- l 21, e7y- l 25, e8y- l 29...
    %%% Deslocamento em z - e1x linha 3, e2x- l 7, e3x- l 11, e4x- l 15,
    %%% e5x- l 19, e6x- l 23, e7x- l 27, e8x- l 31...
    
    % Update waitbar e mensagem
	waitbar(cntwb/steps,f,sprintf('Passo %d de %d - %4.2f %%',cntwb,steps,(cntwb/steps)*100));
    cntwb= cntwb+1;
    %--------------------------
    
end
delete(f)

% Coluna de cada modo
colmodo1= 40;               % indice da freq 205
colmodo2= 99;               % indice da freq 500
colmodo3= 197;              % indice da freq 990
colmodo4= 279;              % indice da freq 1400

%% Resultados
fname = fullfile(pwd, 'output_figures');                                        % Diretoria para salvar as figuras
savefig= 0;                                                                     % savefig=1-> salva; savefig=0 -> nao salva
if savefig && ~exist(fname, 'dir'), mkdir(fname); end

%% Locus do mancal

figure
th_c = linspace(0, 2*pi, 200);                     % reference circles (viscircles is MATLAB-only)
for r_c = [8 6 4 2]*1e-5
    plot(r_c*cos(th_c), r_c*sin(th_c), 'Color', [0.8 0.8 0.8], 'LineWidth', 1); hold on;
end
plot(ey(1,:),ez(1,:),'ko-',ey(2,:),ez(2,:),'bo-','linewidth', 1.5,'MarkerFaceColor','white')
xlim([-0.00006 0.00006])
ylim([-0.00006 0.00006])
set(gcf,'position',[250,150,500,500])
grid on; axis equal
set(gca,'FontSize',14)
legend('Bearing 1','Bearing 2','FontSize',10)
title('Bearing locus','FontWeight','bold')
xlabel('E*cos(\theta)','FontWeight','bold')
ylabel('E*sin(\theta)','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('LocusMancal_caso%d', caso)), 'png');
end

%% Fatores de amortecimento e rigidez dos mancais

% Fat. de rigidez x vel. angular do eixo Mancal 1
figure()
plot(omghz,k11(1,:),omghz,k12(1,:),omghz,k21(1,:),omghz,k22(1,:),'linewidth', 1.5,'MarkerFaceColor','white')
legend('k_1_1','k_1_2','k_2_1','k_2_2','location','southwest','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft speed \Omega [Hz]','FontWeight','bold')
ylabel('k_i_k','FontWeight','bold')
title('Stiffness coefficients k_{ij} vs shaft speed \Omega - bearing 1','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FatKManc1_caso%d', caso)), 'png');
end

% Fat. de rigidez x vel. angular do eixo Mancal 2
figure()
plot(omgrd,k11(2,:),omgrd,k12(2,:),omgrd,k21(2,:),omgrd,k22(2,:),'linewidth', 1.5,'MarkerFaceColor','white')
legend('k_1_1','k_1_2','k_2_1','k_2_2','location','southwest','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft speed \Omega [Hz]','FontWeight','bold')
ylabel('k_i_k','FontWeight','bold')
title('Stiffness coefficients k_{ij} vs shaft speed \Omega - bearing 2','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FatKManc2_caso%d', caso)), 'png');
end

% Fat. de amortecimento x vel. angular do eixo Mancal 1
figure()
plot(omgrd,d11(1,:),omgrd,d12(1,:),omgrd,d21(1,:),'-.',omgrd,d22(1,:),'linewidth', 1.5,'MarkerFaceColor','white')
legend('d_1_1','d_1_2','d_2_1','d_2_2','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft speed \Omega [Hz]','FontWeight','bold')
ylabel('d_i_k','FontWeight','bold')
title('Damping coefficients d_{ij} vs shaft speed \Omega','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FatDManc1_caso%d', caso)), 'png');
end

% Fat. de amortecimento x vel. angular do eixo Mancal 2
figure()
plot(omgrd,d11(2,:),omgrd,d12(2,:),omgrd,d21(2,:),'-.',omgrd,d22(2,:),'linewidth', 1.5,'MarkerFaceColor','white')
legend('d_1_1','d_1_2','d_2_1','d_2_2','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft speed \Omega [Hz]','FontWeight','bold')
ylabel('d_i_k','FontWeight','bold')
title('Damping coefficients d_{ij} vs shaft speed \Omega','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FatDManc2_caso%d', caso)), 'png');
end

%% FRF Mancal 1
figure
plot(omgrd,abs(X(incel(LMc1,1),:)),'.-b',omgrd,abs(X(incel(LMc1,3),:)),'.-k','linewidth', 1,'MarkerFaceColor','white')
legend('y','z','FontSize',10)
xlim([10 1500])
grid on
set(gca,'FontSize',14)
xlabel('Shaft speed [rad/s]','FontWeight','bold')
ylabel('Amplitude','FontWeight','bold')
title('FRF - bearing 1','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FRFMancal1_caso%d', caso)), 'png'); 
end

% FRF Mancal 2
figure
plot(omgrd,abs(X(incel(LMc2,1),:)),'.-b',omgrd,abs(X(incel(LMc2,3),:)),'.-k','linewidth', 1,'MarkerFaceColor','white')
legend('y','z','FontSize',10)
xlim([10 1500])
grid on
set(gca,'FontSize',14)
xlabel('Shaft speed [rad/s]','FontWeight','bold')
ylabel('Amplitude','FontWeight','bold')
title('FRF - bearing 2','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FRFMancal2_caso%d', caso)), 'png'); 
end

% FRF Disco
figure
plot(omgrd,abs(X(incel(LDc,1),:)),'.-b',omgrd,abs(X(incel(LDc,3),:)),'.-k','linewidth', 1,'MarkerFaceColor','white')
legend('y','z','FontSize',10)
xlim([10 1500])
grid on
set(gca,'FontSize',14)
xlabel('Shaft speed [rad/s]','FontWeight','bold')
ylabel('Amplitude','FontWeight','bold')
title('FRF - disk','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FRFdisco_caso%d', caso)), 'png'); 
end

%% Modos 2D
%% Modo 1
figure
plot(cordf,X(1:4:end-3,colmodo1),'o-b',cordf,X(3:4:end-1,colmodo1),'o-k','linewidth', 1.5,'MarkerFaceColor','white'); hold on; grid minor;
plot(cordf,zeros(length(cordf),1), '--k')
legend('y','z','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft position [m]','FontWeight','bold')
ylabel('Amplitude','FontWeight','bold')
title('Deflected shape - 1st mode','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal1_caso%d', caso)), 'png'); 
end

% Modo 2
figure
plot(cordf,X(1:4:end-3,colmodo2),'o-b',cordf,X(3:4:end-1,colmodo2),'o-k','linewidth', 1.5,'MarkerFaceColor','white'); hold on; grid minor;
plot(cordf,zeros(length(cordf),1), '--k')
legend('y','z','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft position [m]','FontWeight','bold')
ylabel('Amplitude','FontWeight','bold')
title('Deflected shape - between 1st and 2nd modes','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal2_caso%d', caso)), 'png'); 
end

% Modo 3
figure
plot(cordf,X(1:4:end-3,colmodo3),'o-b',cordf,X(3:4:end-1,colmodo3),'o-k','linewidth', 1.5,'MarkerFaceColor','white'); hold on; grid minor;
plot(cordf,zeros(length(cordf),1), '--k')
legend('y','z','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft position [m]','FontWeight','bold')
ylabel('Amplitude','FontWeight','bold')
title('Deflected shape - 2nd mode','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal3_caso%d', caso)), 'png'); 
end

% Modo 4
figure
plot(cordf,X(1:4:end-3,colmodo4),'o-b',cordf,X(3:4:end-1,colmodo4),'o-k','linewidth', 1.5,'MarkerFaceColor','white'); hold on; grid minor;
plot(cordf,zeros(length(cordf),1), '--k')
legend('y','z','FontSize',10)
grid on
set(gca,'FontSize',14)
xlabel('Shaft position [m]','FontWeight','bold')
ylabel('Amplitude','FontWeight','bold')
title('Deflected shape - above the 2nd mode','FontWeight','bold')
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal4_caso%d', caso)), 'png'); 
end

%% Modos 3D
%% Modo 1
veczro1= zeros(length(X(1:4:end-3,colmodo1)),1);
Yelps= real(X(1:4:end-3,colmodo1));
Zelps= real(X(3:4:end-1,colmodo1));

for i=1:length(cordf)               % Interpolando as elipses de cada no do eixo
    [Y(i,:),Z(i,:)]= intrpElipse(Yelps(i),Zelps(i));
    cordfelps(i,1:length(Y(i,:)))= cordf(i);
end

figure();
plot3(cordf,veczro1,X(1:4:end-3,colmodo1), '--r', 'linewidth', 1.4); hold on; grid on;
plot3(cordf,X(3:4:end-1,colmodo1),veczro1, '--b', 'linewidth', 1.4)
plot3(cordf,veczro1,veczro1, '--k', 'linewidth', 0.8)
for i=1:length(cordf)                                               
    plot3(cordfelps(i,:),Z(i,:),Y(i,:), 'k', 'linewidth', 0.2)                  % Plotando as elipses de cada no do eixo
end
grid on; grid minor;
set(gcf,'position',[250,150,1300,800])
set(gca,'FontSize',18)
legend('Displacement y','Displacement z','Shaft axis','Node whirl orbit','FontSize',14,'location','northeast')
title('Mode 1 - 3D whirl','FontWeight','bold')
xlabel('x','FontWeight','bold')
ylabel('z','FontWeight','bold')
zlabel('y','FontWeight','bold')
view(-20,20)
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal3D1_caso%d', caso)), 'png'); 
   saveas(gca, fullfile(fname, sprintf('FModal3D1_caso%d', caso)), 'fig'); 
end

% Modo 2
veczro1= zeros(length(X(1:4:end-3,colmodo2)),1);
Yelps= real(X(1:4:end-3,colmodo2));
Zelps= real(X(3:4:end-1,colmodo2));

for i=1:length(cordf)               % Interpolando as elipses de cada no do eixo
    [Y(i,:),Z(i,:)]= intrpElipse(Yelps(i),Zelps(i));    
    cordfelps(i,1:length(Y(i,:)))= cordf(i);
end

figure();
plot3(cordf,veczro1,X(1:4:end-3,colmodo2), '--r', 'linewidth', 1.4); hold on; grid on;
plot3(cordf,X(3:4:end-1,colmodo2),veczro1, '--b', 'linewidth', 1.4)
plot3(cordf,veczro1,veczro1, '--k', 'linewidth', 0.8)
for i=1:length(cordf)                                               
    plot3(cordfelps(i,:),Z(i,:),Y(i,:), 'k', 'linewidth', 0.2)                  % Plotando as elipses de cada no do eixo
end
grid on; grid minor;
set(gcf,'position',[250,150,1300,800])
set(gca,'FontSize',18)
legend('Displacement y','Displacement z','Shaft axis','Node whirl orbit','FontSize',14,'location','northeast')
title('Mode 2 - 3D whirl','FontWeight','bold')
xlabel('x','FontWeight','bold')
ylabel('z','FontWeight','bold')
zlabel('y','FontWeight','bold')
view(-20,20)
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal3D2_caso%d', caso)), 'png'); 
   saveas(gca, fullfile(fname, sprintf('FModal3D2_caso%d', caso)), 'fig'); 
end

% Modo 3
veczro1= zeros(length(X(1:4:end-3,colmodo3)),1);
Yelps= real(X(1:4:end-3,colmodo3));
Zelps= real(X(3:4:end-1,colmodo3));

for i=1:length(cordf)               % Interpolando as elipses de cada no do eixo
    [Y(i,:),Z(i,:)]= intrpElipse(Yelps(i),Zelps(i));    
    cordfelps(i,1:length(Y(i,:)))= cordf(i);
end

figure();
plot3(cordf,veczro1,X(1:4:end-3,colmodo3), '--r', 'linewidth', 1.4); hold on; grid on;
plot3(cordf,X(3:4:end-1,colmodo3),veczro1, '--b', 'linewidth', 1.4)
plot3(cordf,veczro1,veczro1, '--k', 'linewidth', 0.8)
for i=1:length(cordf)                                               
    plot3(cordfelps(i,:),Z(i,:),Y(i,:), 'k', 'linewidth', 0.2)                  % Plotando as elipses de cada no do eixo
end
grid on; grid minor;
set(gcf,'position',[250,150,1300,800])
set(gca,'FontSize',18)
legend('Displacement y','Displacement z','Shaft axis','Node whirl orbit','FontSize',14,'location','northeast')
title('Mode 3 - 3D whirl','FontWeight','bold')
xlabel('x','FontWeight','bold')
ylabel('z','FontWeight','bold')
zlabel('y','FontWeight','bold')
view(-20,20)
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal3D3_caso%d', caso)), 'png'); 
   saveas(gca, fullfile(fname, sprintf('FModal3D3_caso%d', caso)), 'fig'); 
end

% Modo 4
veczro1= zeros(length(X(1:4:end-3,colmodo4)),1);
Yelps= real(X(1:4:end-3,colmodo4));
Zelps= real(X(3:4:end-1,colmodo4));

for i=1:length(cordf)               % Interpolando as elipses de cada no do eixo
    [Y(i,:),Z(i,:)]= intrpElipse(Yelps(i),Zelps(i));    
    cordfelps(i,1:length(Y(i,:)))= cordf(i);
end

figure();
plot3(cordf,veczro1,X(1:4:end-3,colmodo4), '--r', 'linewidth', 1.4); hold on; grid on;
plot3(cordf,X(3:4:end-1,colmodo4),veczro1, '--b', 'linewidth', 1.4)
plot3(cordf,veczro1,veczro1, '--k', 'linewidth', 0.8)
for i=1:length(cordf)                                               
    plot3(cordfelps(i,:),Z(i,:),Y(i,:), 'k', 'linewidth', 0.2)                  % Plotando as elipses de cada no do eixo
end
grid on; grid minor;
set(gcf,'position',[250,150,1300,800])
set(gca,'FontSize',18)
legend('Displacement y','Displacement z','Shaft axis','Node whirl orbit','FontSize',14,'location','northeast')
title('Mode 4 - 3D whirl','FontWeight','bold')
xlabel('x','FontWeight','bold')
ylabel('z','FontWeight','bold')
zlabel('y','FontWeight','bold')
view(-20,20)
if savefig==1
   saveas(gca, fullfile(fname, sprintf('FModal3D4_caso%d', caso)), 'png'); 
   saveas(gca, fullfile(fname, sprintf('FModal3D4_caso%d', caso)), 'fig'); 
end

