% Analise de maquinas rotativas - Exercicio 4 - Mancais hidrodinamicos
% Nicolas Londero Gass
clear all; close all; clc
warning off;

%% Dados iniciais do mancal

dm= 30*10^-3;                                       % Diametro do mancal
cm= 20*10^-3;                                       % Comprimento mancal
frd= 0.09*10^-3;                                    % Folga Radial
eta= 7*10^-2;                                       % Viscosidade absoluta
R= dm/2;                                            % Raio do mancal

omg= 9.5*(2*pi);                                    % Omega, tirado do slide 39, aula 10\11
E= 0.1;                                             % Excentricidade, tirado do slide 39, aula 10\11
Sm= (pi/2)*(E/(1-E^2)^2)*sqrt(1-E^2+((4/pi)*E)^2);  % S* - Sommerfeld modificado
Fn= (eta*cm^3*omg*R)/(2*frd^2);                     % Força de fricção
F0= Fn*Sm;                                          % Carga estática aplicada no mancal

omghz= 0.2:0.1:50;                                  % Vetor omega em Hertz
omgrd= omghz*(2*pi);                                % Vetor omega em radianos

valest= [0.6; 0.1];                                 % Valores estimados de epsilon e phi, tirado do slide 39

n= size(omgrd,2);

% Inicializando vetores
Phi= zeros(1,n); phi= zeros(1,n); Sm22= zeros(1,n); 
Sm2= zeros(1,n); g22= zeros(1,n); g21= zeros(1,n); 
g12= zeros(1,n); g11= zeros(1,n); b22= zeros(1,n); 
b21= zeros(1,n); b12= zeros(1,n); b11= zeros(1,n); 
k11= zeros(1,n); k12= zeros(1,n); k21= zeros(1,n); 
k22= zeros(1,n); d11= zeros(1,n); d12= zeros(1,n); 
d21= zeros(1,n); d22= zeros(1,n); E3= zeros(1,n); 

% Waitbar ------------------------------------------
steps = length(omgrd);
f = waitbar(0,sprintf('Passo %d de %d - %5.2f %%',0,steps,(0)*100));
setappdata(f,'canceling',0);
% --------------------------------------------------

tic
for i=1:length(omgrd)
    
    if getappdata(f,'canceling')
        break
    end
    
    K0= (eta*omgrd(i)*R*(cm^3))/(frd^2);        % Fator de forca do mancal curto
    [E3(i),Phi(i)]= MNR_bearing(K0, F0, valest);      % Epsilon e Phi por Newton-Raphson (Jacobiano analitico)
    
    valest= [E3(i); Phi(i)];                                                % Variaveis finais
    phi(i)= Phi(i)*(180/pi);                                                % phi em graus
    
    Fn2= (eta*cm^3*omgrd(i)*(R))/(2*frd^2);                                 % Força de fricção para omghz(i)
    
    Sm2(i)= F0/Fn2;                                                         % Sommerfeld modificado
    Sm22(i)= (pi/2)*(E3(i)/(1-E3(i)^2)^2)*sqrt(1-E3(i)^2+((4/pi)*E3(i))^2); % Sommerfeld modificado
    
    % Calculos de Kik e Cik e variaveis auxiliares
    
    A= 4/(pi^2+(16-pi^2)*E3(i)^2)^(3/2);                                      
    
    g11(i)= (2*pi^2+(16-pi^2)*E3(i)^2)*A;                                                           % Gamma11
    g12(i)= (pi/4)* ((pi^2-(2*pi^2)*E3(i)^2-(16-pi^2)*E3(i)^4)/ (E3(i)*sqrt(1-E3(i)^2)))*A;         % Gamma12
    g21(i)= (-pi/4)* ((pi^2+(32+pi^2)*(E3(i)^2)+(32-2*pi^2)*(E3(i)^4))/ (E3(i)*sqrt(1-E3(i)^2)))*A; % Gamma21
    g22(i)= (((pi^2)+(32+(pi^2))*(E3(i)^2)+ (32-2*(pi^2))*(E3(i)^4))/ (1-(E3(i)^2)))*A;              % Gamma22 (denominador 1-eps^2)
    
    b11(i)= (pi/2)*((sqrt(1-E3(i)^2))/E3(i))* (pi^2+(2*pi^2-16)*E3(i)^2)*A;                         % Beta11
    b12(i)= -(2*pi^2+(4*pi^2-32)*E3(i)^2)*A;                                                        % Beta12
    b21(i)= b12(i);                                                                                 % Beta21
    b22(i)= (pi/2)*((pi^2+(48-2*pi^2)*E3(i)^2+(pi^2)*E3(i)^4)/ (E3(i)*sqrt(1-E3(i)^2)))*A;          % Beta22
    
    k11(i)= g11(i)*(F0/frd);                        % k11
    k12(i)= g12(i)*(F0/frd);                        % k12
    k21(i)= g21(i)*(F0/frd);                        % k21
    k22(i)= g22(i)*(F0/frd);                        % k22
    
    d11(i)= b11(i)*(F0/(frd*omgrd(i)));             % d11
    d12(i)= b12(i)*(F0/(frd*omgrd(i)));             % d12
    d21(i)= b21(i)*(F0/(frd*omgrd(i)));             % d21
    d22(i)= b22(i)*(F0/(frd*omgrd(i)));             % d22
    
    % Update waitbar e mensagem
    waitbar(i/steps,f,sprintf('Passo %d de %d - %4.2f %%',i,steps,(i/steps)*100));
    
end
tt(1)= toc;

delete(f)

kij= [k11;k12;k21;k22];
dij= [d11;d12;d21;d22];

%% Rotor Laval com mancais lubrificados

De= 10*10^(-3);                                     % Diametro do eixo [m]
Ce= 800*10^(-3);                                    % Comprimento do eixo [m]
I= (pi*(De/2)^4)/4;                                 % Momento de inercia eixo [m^4]

rho= 7850;                                          % Densidade do aco [kg/m³]
Dd= 90*10^(-3);                                     % Diametro do disco [m]
Cd= 20*10^(-3);                                     % Comprimento do disco [m]
Ad= pi*(Dd/2)^2;                                    % Area do disco [m²]
Vd= Ad*Cd;                                          % Volume do disco [m³]
md= Vd*rho;                                         % Massa do disco [kg]
% Jg= (md*(Dd^2))/8;                                  % Momento de inercia polar [m^4]

ME= 200*10^9;                                       % Modulo de elasticidade do aco [Pa]
ps= 0.3;                                            % Coef Poisson [kg/m³]
k= (48*ME*I)/(Ce^3);                                % Constante elastica do eixo [N/m]
c= 10^(-4)*k;                                       % Coef de amortecimento [N.s/m]
me= 10^(-4);                                        % Desbalanceamento [kg.m]
g= 9.81;                                            % Gravidade [m/s²]

cc= [0 0 0 0 0 0];                                  % Condicoes de contorno
vars= [k;c;g;md;me];

% TF= 10;
% dt= TF/steps;
tvec= 0:0.001:10;                                   % Vetor de tempo
N= length(tvec);                                    % Tamanho do vetor de tempo

uy= zeros(length(omgrd),N);
uz= zeros(length(omgrd),N);
ym= zeros(length(omgrd),N);
zm= zeros(length(omgrd),N);

% Waitbar ------------------------------------------
steps = length(omgrd);
f = waitbar(0,sprintf('Passo %d de %d - %5.2f %%',0,steps,(0)*100));
setappdata(f,'canceling',0);
% --------------------------------------------------

tic
for i=1:length(omgrd)
    
    if getappdata(f,'canceling')
        break
    end
    
    [t, res]= ode45(@(t,x)eqsd_v1(t,x,omgrd(i),kij(:,i),dij(:,i),vars), tvec, cc);
    
    uy(i,:)= res(:,1);                              % Separando variaveis
    uz(i,:)= res(:,2);                              %
    ym(i,:)= res(:,3);                              %
    zm(i,:)= res(:,4);                              %
    
    muy= mean(uy(i,round(N*0.98):end));             % Media dos valores
    muz= mean(uz(i,round(N*0.98):end));             %
    mym= mean(ym(i,round(N*0.98):end));             %
    mzm= mean(zm(i,round(N*0.98):end));             %
    
    uy2= uy(i,round(N*0.98):end)-muy;               %
    uz2= uz(i,round(N*0.98):end)-muz;               %
    ym2= ym(i,round(N*0.98):end)-mym;               %
    zm2= zm(i,round(N*0.98):end)-mzm;               %
    
    xd(i)= max(sqrt(uy2.^2+uz2.^2));                % Deslocamento do disco
    xm(i)= max(sqrt(ym2.^2+zm2.^2));                % Deslocamento do mancal
    
    % Update waitbar e mensagem
    waitbar(i/steps,f,sprintf('Passo %d de %d - %4.2f %%',i,steps,(i/steps)*100));
    
end
tt(2)=toc;

delete(f)

ey= (E3*frd/2).*sin(phi/(180/pi));
ez= -(E3*frd/2).*cos(phi/(180/pi));

fprintf('Tempo de calculo de K e C: %.2fs \nTempo de calculo dos deslocamentos: %.2fs \n', tt(1), tt(2));

%% Resultados

save=0;                                             % save=1 -> Salva as figuras
if save == 1
    % Definindo a pasta para salvar as figuras em PNG
    fname = fullfile(pwd, 'output_figures');  % Output folder for figures
    if ~exist(fname, 'dir'), mkdir(fname); end
end

% Locus (centro do eixo em relação à folga do mancal)
figure(1)
if exist('polarplot') ~= 0                          % MATLAB (R2016a+)
    polarplot(phi*(pi/180),E3,'linewidth', 1.5); hold on
    polarplot(phi(1)*(pi/180),E3(1),'rx','linewidth', 1.5, 'markersize', 9)
    polarplot(phi(end)*(pi/180),E3(end),'kx','linewidth', 1.5, 'markersize', 9)
    grid minor
    set(gca,'FontSize',12)
    thetalim([0 90])
    rlim([0 1])
    rticks([0 0.2 0.4 0.6 0.8 1])
    thetaticks([0 45 90])
    rticklabels({'\epsilon','0.2','0.4','0.6','0.8','1'})
else                                                % Octave: same locus in Cartesian axes
    th= phi*(pi/180);
    plot(E3.*cos(th), E3.*sin(th),'linewidth', 1.5); hold on
    plot(E3(1)*cos(th(1)), E3(1)*sin(th(1)),'rx','linewidth', 1.5, 'markersize', 9)
    plot(E3(end)*cos(th(end)), E3(end)*sin(th(end)),'kx','linewidth', 1.5, 'markersize', 9)
    tc= linspace(0, pi/2, 100);
    for rr= 0.2:0.2:1                               % reference arcs (epsilon = 0.2 ... 1)
        plot(rr*cos(tc), rr*sin(tc), ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off')
    end
    axis equal; xlim([0 1]); ylim([0 1]); grid on
    set(gca,'FontSize',12)
    xlabel('\epsilon cos\phi'); ylabel('\epsilon sin\phi')
end
title({'Bearing locus'})
legend('Locus','Start point','End point')
if save == 1
    saveas(gca, fullfile(fname, 'LocusMancal'), 'png');
end

% Excentricidade x vel. angular
figure(2)
plot(omghz,E3,'linewidth', 1.5)
ylim([0 0.71])
grid on
set(gca,'FontSize',12)
legend('\epsilon')
title(['Eccentricity \epsilon vs shaft speed \Omega'],'FontWeight','bold')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Eccentricity \epsilon')
if save == 1
    saveas(gca, fullfile(fname, 'excentricidadeXVelAng'), 'png');
end

% Sommerfeld e phi x excentricidade
figure(3)
[ax3, hS, hA] = plotyy(E3, Sm2, E3, phi, 'semilogy', 'plot');   % plotyy: MATLAB and Octave
set([hS hA], 'linewidth', 1.5)
xlim(ax3(1), [0 0.8]); xlim(ax3(2), [0 0.8])
grid(ax3(1), 'on')
set(ax3, 'FontSize', 12)
xlabel(ax3(1), 'Eccentricity \epsilon')
ylabel(ax3(1), 'Modified Sommerfeld number S*')
ylabel(ax3(2), '\alpha [deg]')
legend([hS hA], 'S*', '\alpha', 'location', 'southwest')
title(['S* and attitude angle \alpha vs eccentricity \epsilon'],'FontWeight','bold')
if save == 1
    saveas(gca, fullfile(fname, 'SommPhiXExcentr'), 'png');
end

% Gamma x excentricidade
figure(4)
semilogy(E3,g11,E3,g12,E3,-g21,E3,g22,'linewidth', 1.5,'MarkerFaceColor','white')
xlim([0 0.7])
grid on
set(gca,'FontSize',12)
xlabel('Eccentricity \epsilon')
ylabel('\gamma_i_k')
title(['Stiffness coefficients \gamma_{ij} vs eccentricity \epsilon'],'FontWeight','bold')
legend('\gamma_{11}','\gamma_{12}','-\gamma_{21}','\gamma_{22}', 'location', 'southwest')
if save == 1
    saveas(gca, fullfile(fname, 'GammaXExcentr'), 'png');
end

% Beta x excentricidade
figure(5)
semilogy(E3,b11,E3,-b12,E3,-b21,'-.',E3,b22,'linewidth', 1.5,'MarkerFaceColor','white')
xlim([0 0.71])
grid on
set(gca,'FontSize',12)
xlabel('Eccentricity \epsilon')
ylabel('\beta_i_k')
title(['Damping coefficients \beta_{ij} vs eccentricity \epsilon'],'FontWeight','bold')
legend('\beta_{11}','-\beta_{12}','-\beta_{21}','\beta_{22}', 'location', 'northeast')
if save == 1
    saveas(gca, fullfile(fname, 'Beta'), 'png');
end

% Fat. de rigidez x vel. angular do eixo
figure(6)
plot(omghz,k11,omghz,k12,omghz,k21,omghz,k22,'linewidth', 1.5,'MarkerFaceColor','white')
legend('k_1_1','k_1_2','k_2_1','k_2_2','location','southwest')
grid on
set(gca,'FontSize',12)
xlabel('Shaft speed \Omega [Hz]')
ylabel('k_i_k')
title(['Stiffness coefficients k_{ij} vs shaft speed \Omega'],'FontWeight','bold')
if save == 1
    saveas(gca, fullfile(fname, 'FatRigidez'), 'png');
end

% Fat. de amortecimento x vel. angular do eixo
figure(7)
plot(omghz,d11,omghz,d12,omghz,d21,'-.',omghz,d22,'linewidth', 1.5,'MarkerFaceColor','white')
legend('d_1_1','d_1_2','d_2_1','d_2_2')
grid on
set(gca,'FontSize',12)
xlabel('Shaft speed \Omega [Hz]')
ylabel('d_i_k')
title(['Damping coefficients d_{ij} vs shaft speed \Omega'],'FontWeight','bold')
if save == 1
    saveas(gca, fullfile(fname, 'FatAmortec'), 'png');
end

% Orbita do disco
figure(8)
plot(uy(49,round(N*0.98):end),uz(49,round(N*0.98):end),uy(149,round(N*0.98):end),uz(149,round(N*0.98):end),uy(249,round(N*0.98):end),uz(249,round(N*0.98):end),...
     'linewidth', 1.5,'MarkerFaceColor','white')
%     uy(305,round(N*0.98):end),uz(305,round(N*0.98):end),'linewidth', 1.5,'MarkerFaceColor','white')
legend('\Omega = 5 Hz','\Omega = 15 Hz','\Omega = 25 Hz')
grid on
ylim([-4e-3 2e-3])
xlim([-4e-3 4e-3])
set(gca,'FontSize',12)
xlabel('Disk displacement y [m]')
ylabel('Disk displacement z [m]')
title(['Disk centre orbit'],'FontWeight','bold')
if save == 1
    saveas(gca, fullfile(fname, 'orbitaDisco'), 'png');
end

% Orbita da folga do mancal
xfm= 0:0.1:2*pi+0.1;
xcirc= (frd/2)*sin(xfm);
ycirc= (frd/2)*cos(xfm);

% Orbita do mancal
figure(9)
hold on
grid on
plot(ym(49,round(N*0.98):end),zm(49,round(N*0.98):end),ym(149,round(N*0.98):end),zm(149,round(N*0.98):end),...
     ym(249,round(N*0.98):end),zm(249,round(N*0.98):end),ey,ez,'k','linewidth', 1.5,'MarkerFaceColor','white')
plot(xcirc,ycirc','linewidth', 1.5,'MarkerFaceColor','white')
ylim([-5e-5 5e-5])
xlim([-6.5e-5 6.5e-5])
set(gca,'FontSize',12)
xlabel('Journal displacement y [m]')
ylabel('Journal displacement z [m]')
title(['Journal centre orbit'],'FontWeight','bold')
legend('\Omega = 5 Hz','\Omega = 15 Hz','\Omega = 25 Hz','Locus','Bearing clearance')
if save == 1
    saveas(gca, fullfile(fname, 'orbitaMancal'), 'png');
end

% FRF do disco
figure(10)
plot(omghz,xd,'linewidth', 1.5,'MarkerFaceColor','white')
grid on
xlim([0 30])
set(gca,'FontSize',12)
legend('FRF')
xlabel('Frequency [Hz]')
ylabel('Displacement [m]')
title(['FRF - disk centre'],'FontWeight','bold')
if save == 1
    saveas(gca, fullfile(fname, 'FrfDisco'), 'png');
end

% FRF do mancal
figure(11)
plot(omghz,xm,'linewidth', 1.5,'MarkerFaceColor','white')
grid on
xlim([0 30])
set(gca,'FontSize',12)
legend('FRF')
xlabel('Frequency [Hz]')
ylabel('Displacement [m]')
title(['FRF - journal centre'],'FontWeight','bold')
if save == 1
    saveas(gca, fullfile(fname, 'FrfMancal'), 'png');
end
