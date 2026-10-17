%% Elementos Finitos 3 - Análise de estabilidade dos métodos de integração - pelo metodo da Matriz A
% Nicolas Londero Gass - nicolas.gass@eac.ufsm.br
clear all; close all; clc;
warning off;

%% Parametros iniciais

omg= 1;
alp= 0.25;
delt= 0.5;
tht1= 1.4;
tht2= 1.37;
tht3= 1.366;
thtv1= 1:0.1:1.3;                                   % Theta variavel
thtv2= 1.3:0.01:1.36;                               % Theta variavel
thtv3= 1.361:0.001:1.366;                           % Theta variavel
thtv4= 1.366:0.00001:1.3662;                        % Theta variavel
thtv5= 1.3662:0.01:1.5;                             % Theta variavel
thtv6= 1.37;                                        % Theta variavel
thtv= [thtv1, thtv2, thtv3, thtv4, thtv5, thtv6];   % Vetor theta final

thtb1= 0.01:0.01:1.36;                              % Theta1 grafico 2
thtb2= 1.366:0.00001:1.368;                         % Theta2 grafico 2
thtb3= 1.37:0.01:3.4;                               % Theta3 grafico 2
thtb= [thtb1, thtb2, thtb3];                        % Theta grafico 2

dtmin= 0.0001;
dtmax= 10000;
dt1= dtmin:0.001:10;
dt2= 10:0.1:1000;
dt3= 1000:10:dtmax;
dtf= [dt1,dt2,dt3];                                 % Vetor dt final (discretizacao variavel)

% Alocando vetores
pd= zeros(1,length(dtf)); ph= zeros(1,length(dtf)); pn= zeros(1,length(dtf));
pw1= zeros(1,length(dtf)); pw2= zeros(1,length(dtf)); pw3= zeros(1,length(dtf));
pwv= zeros(length(thtv),length(dtf)); pwv2= zeros(1,length(thtb));

for i=1:length(dtf)
    % 
    dt= dtf(i)*(2*pi)/omg;
    
    % Matriz A - Diferenca central V
    Ad= [(2-omg^2*dt^2)/1, -1;
                1,          0];

    % Matriz A - Houbolt V
    beth= (2/((omg^2)*(dt^2))+1)^(-1);
    
    Ah= [(5*beth)/((omg^2)*(dt^2)),   -((4*beth)/((omg^2)*(dt^2))),   ((beth)/((omg^2)*(dt^2)));
                   1,                               0,                            0;
                   0,                               1,                            0];
              
    % Matriz A - Wilson theta (tht=1.4)
    betw1= ((tht1/((omg*dt)^2))+(tht1^3)/(6))^(-1);
    
    Aw1= [          1-(betw1*tht1^2)/3-(1/tht1),                (1/dt)*(-betw1*tht1),      (1/dt^2)*(-betw1);
                dt*(1-(1/(2*tht1))-((betw1*tht1^2)/6)),          1-((betw1*tht1)/2),       (1/dt)*(-betw1/2);
           (dt^2)*((1/2)-(1/(6*tht1))-((betw1*tht1^2)/18)),    dt*(1-((betw1*tht1)/6)),       1-(betw1/6)];

	% Matriz A - Wilson theta (tht=1.37)
    betw2= ((tht2/((omg^2)*(dt^2)))+(tht2^3)/(6))^(-1);
    
    Aw2= [           1-(betw2*tht2^2)/3-1/tht2,                (1/dt)*(-betw2*tht2),         (1/dt^2)*(-betw2);
               dt*(1-(1/(2*tht2))-((betw2*tht2^2)/6)),          1-((betw2*tht2)/2),          (1/dt)*(-betw2/2);
          (dt^2)*((1/2)-(1/(6*tht2))-((betw2*tht2^2)/18)),    dt*(1-((betw2*tht2)/6)),         (1)-(betw2/6)];

	% Matriz A - Wilson theta (tht=1.366)
    betw3= ((tht3/((omg^2)*(dt^2)))+(tht3^3)/(6))^(-1);
    
    Aw3= [           1-(betw3*tht3^2)/3-1/tht3,                (1/dt)*(-betw3*tht3),      (1/dt^2)*(-betw3);
               dt*(1-(1/(2*tht3))-((betw3*tht3^2)/6)),          1-((betw3*tht3)/2),       (1/dt)*(-betw3/2);
          (dt^2)*((1/2)-(1/(6*tht3))-((betw3*tht3^2)/18)),    dt*(1-((betw3*tht3)/6)),      (1)-(betw3/6)];

    % Matriz A - Newmark
    betn= ((1/((omg^2)*(dt^2)))+alp)^(-1);
    
    An= [         -((1/2)-alp)*betn,               (1/dt)*(-betn),      (1/dt^2)*(-betn);
            dt*(1-delt-(1/2-alp)*delt*betn),         1-betn*delt,      (1/dt)*(-betn*delt);
         (dt^2)*((1/2)-alp-(1/2-alp)*alp*betn),    dt*(1-betn*alp),       (1-betn*alp)];
     
    % Maximo do absoluto dos autovalores das matrizes A
    pd(i)= max(abs(eig(Ad)));
    ph(i)= max(abs(eig(Ah)));
    pn(i)= max(abs(eig(An)));
    pw1(i)= max(abs(eig(Aw1)));
    pw2(i)= max(abs(eig(Aw2)));
    pw3(i)= max(abs(eig(Aw3)));
        
	for j= 1:length(thtv)

        % Matriz A - Wilson theta (tht variavel - pts de interesse)
        betwv= ((thtv(j)/((omg*dt)^2))+(thtv(j)^3)/(6))^(-1);

        Awv= [          1-(betwv*thtv(j)^2)/3-(1/thtv(j)),                (1/dt)*(-betwv*thtv(j)),           (1/dt^2)*(-betwv);
                   dt*(1-(1/(2*thtv(j)))-((betwv*thtv(j)^2)/6)),          1-((betwv*thtv(j))/2),            (1/dt)*(-betwv/2);
               (dt^2)*((1/2)-(1/(6*thtv(j)))-((betwv*thtv(j)^2)/18)),    dt*(1-((betwv*thtv(j))/6)),            1-(betwv/6)];

        pwv(j,i)= max(abs(eig(Awv)));
    end
    
    if i==length(dtf)
        dt=10000;
        for j=1:length(thtb)
            % Matriz A - Wilson theta (tht variavel)
            betwv= ((thtb(j)/((omg*dt)^2))+(thtb(j)^3)/(6))^(-1);

            Awv= [          1-(betwv*thtb(j)^2)/3-(1/thtb(j)),                (1/dt)*(-betwv*thtb(j)),           (1/dt^2)*(-betwv);
                       dt*(1-(1/(2*thtb(j)))-((betwv*thtb(j)^2)/6)),          1-((betwv*thtb(j))/2),            (1/dt)*(-betwv/2);
                   (dt^2)*((1/2)-(1/(6*thtb(j)))-((betwv*thtb(j)^2)/18)),    dt*(1-((betwv*thtb(j))/6)),            1-(betwv/6)];

            pwv2(j)= max(abs(eig(Awv)));
%             disp('if test \n')
        end
    end
end

%% Resultado dos metodos

figure(1)
semilogx(dtf, pd, 'linewidth', 1.6); hold on; grid minor;
semilogx(dtf, ph, 'linewidth', 1.6);
semilogx(dtf, pw1, 'linewidth', 1.6);
semilogx(dtf, pw2, 'linewidth', 1.6);
semilogx(dtf, pn, 'linewidth', 1.6);
axis([dtmin dtmax -0.015 1.4]);
legend('\rho central diff. (+\infty)','\rho Houbolt','\rho Wilson \theta=1.4','\rho Wilson \theta=1.37',...
       '\rho Newmark (\delta=0.5 \alpha=0.25)', 'location', 'southwest', 'fontsize', 11);
xlabel('\Deltat/T', 'fontsize', 13); ylabel('\rho(A)', 'fontsize', 13);
title('Spectral radius of the integration methods (\zeta = 0)', 'fontsize', 15)

%% Grafico 2 Bathe

figure(2);
fig2= plot(thtb, pwv2, 'linewidth', 1.6); hold on; grid on;
plot(thtb(139), pwv2(139), 'rx', 'linewidth', 1.2, 'markersize', 12)
plot([0 3.4], [1 1], 'k', 'linewidth', 1.2);       % yline (portable)
axis([0 3.4 -0.1 2]);
legend('\rho(\theta)','\rho=1', 'location', 'northeast', 'fontsize', 11);
xlabel('\theta', 'fontsize', 13); ylabel('\rho(\theta)', 'fontsize', 13);
title('Spectral radius - Wilson \theta (0<\theta<3.4, \Deltat/T = 10000)', 'fontsize', 15);

if exist('datatip') ~= 0, dattp= datatip(fig2,1.36602,1.00003); end   % MATLAB only

%% Variacao de theta

figure(3)
semilogx(dtf, pwv(1,:), 'linewidth', 1.6); hold on; grid on;
semilogx(dtf, pwv(3,:), 'linewidth', 1.6);
semilogx(dtf, pwv(5,:), 'linewidth', 1.6);
semilogx(dtf, pwv(20,:), 'linewidth', 1.6);
semilogx(dtf, pwv(43,:), 'linewidth', 1.6);
semilogx(dtf, pwv(52,:), 'linewidth', 1.6);
axis([dtmin dtmax 0.5 2.01]);
legend('\theta = 1', '\theta = 1.2', '\theta = 1.3', '\theta = 1.36602', '\theta = 1.4',...
       '\theta = 1.5', 'location', 'northwest', 'fontsize', 11);
xlabel('\Deltat/T', 'fontsize', 13); ylabel('\rho(A)', 'fontsize', 13);
title('Spectral radius - Wilson \theta (1<\theta<1.5)', 'fontsize', 15)
