% Analise de maquinas rotativas - Exercicio 2
% Nicolas Londero Gass
clear all; close all; clc

%% Caso disco descentralizado com mancais rigidos
%% Parametros iniciais
E= 200*(10^9);                          %Modulo de elasticidade [Pa]
rho= 7850;                              %Densidade [kg/m³]
Ps= 0.3;                                %Poisson

%Eixo
De= 0.01;                               %Diametro eixo
Ce= 0.80;                               %Comprimento eixo
a= (2/3)*Ce;                            %Distancia disco-mancal A
b= (1/3)*Ce;                            %Distancia disco-mancal B
Ie= (pi*(De/2)^4)/4;                    %Momento de flexibilidade do eixo [m^4]

%Disco
Dd= 0.09;                               %Diametro disco
Cd= 0.02;                               %Comprimento disco
Ad= pi*(Dd/2)^2;                        %Area do disco [m²]
Vd= Ad*Cd;                              %Volume do disco [m³]
md= Vd*rho;                             %Massa do disco []
Ip= (md*(Dd^2))/8;                      %Momento de inercia polar do disco [m^4]
Id= Ip/2 + (md*Cd^2)/12;                %Momento de inercia diametral do disco [m^4]

Ke= (3*E*Ie*(a^3 + b^3))/(a^3 * b^3);	%Rigidez do eixo (k22)
k= (3*E*Ie*Ce*(a-b))/(a^2 * b^2);       %Rigidez do eixo (k23)
Kr= (3*E*Ie*Ce)/(a*b);                  %Rigidez do eixo (k33)

M= [md 0 0 0;...
    0 md 0 0;...
    0 0 Id 0;...
    0 0 0 Id];

K= [Ke 0 0 -k;...
    0 Ke k 0;...
    0 k Kr 0;...
    -k 0 0 Kr];

G= [0 0 0 0;...
    0 0 0 0;...
    0 0 0 Ip;...
    0 0 -Ip 0];

Z= zeros(size(M));

%%
omg= 0:pi:500*2*pi;
omghz= omg/2/pi;

frq1=[]; frq2=[];
frq3=[]; frq4=[];
for i= 1:length(omg)
    
    GM= [G*omg(i) M;...
          M Z];
    
    KM= [K Z;...
         Z -M];
        
    res= eig(KM,-GM);
    omgn= sort(imag(res))/2/pi;
    
    frq1= [frq1 omgn(5)];
    frq2= [frq2 omgn(6)];
    frq3= [frq3 omgn(7)];
    frq4= [frq4 omgn(8)];
end

%Gráficos caso 1
figure(1)
plot(omghz, omghz,'b--',omghz, omghz*2,'r--','linewidth', 1.5); hold on
plot(omghz, frq1,'m',omghz, frq2,'k',omghz, frq3,'g',omghz, frq4,'b','linewidth', 1.5)
grid on; grid minor;
set(gca,'fontsize',11)
legend({'\omega=\Omega','\omega=2\Omega','Mode 1 (BW)','Mode 1 (FW)','Mode 2 (BW)','Mode 2 (FW)'},'location','northwest')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Frequency \omega [Hz]')
title('Campbell diagram - off-centre disk, rigid bearings')

figure(2)
plot(omghz, omghz,'b--',omghz, omghz*2,'r--','linewidth', 1.5); hold on
plot(omghz, frq1,'m',omghz, frq2,'k',omghz, frq3,'g',omghz, frq4,'b','linewidth', 1.5)
grid on; grid minor; 
ylim([0 50]);
set(gca,'fontsize',11)
legend({'\omega=\Omega','\omega=2\Omega','Mode 1 (BW)','Mode 1 (FW)'},'location','northeast')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Frequency \omega [Hz]')
title('Campbell diagram - off-centre disk, rigid bearings')

figure(3)
plot(omghz, omghz,'b--',omghz, omghz*2,'r--','linewidth', 1.5); hold on
plot(omghz, frq1,'m',omghz, frq2,'k',omghz, frq3,'g',omghz, frq4,'b','linewidth', 1.5)
grid on; grid minor; 
ylim([0 500]);
xlim([0 200]);
set(gca,'fontsize',11)
legend({'\omega=\Omega','\omega=2\Omega','Mode 1 (BW)','Mode 1 (FW)'},'location','northeast')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Frequency \omega [Hz]')
title('Campbell diagram - off-centre disk, rigid bearings')



%% Caso disco descentralizado com mancal flexivel
%% Parametros iniciais
Kmy= Ke;                                %Coeficiente de rigidez do mancal em y
Kmz= 10*Ke;                             %Coeficiente de rigidez do mancal em z
alp= a/Ce;                              %Coeficiente alpha
bet= b/Ce;                              %Coeficiente beta

%Coeficientes de influencia, mancais rigidos e eixo flexivel
h22= (a^2*b^2)/(3*E*Ie*Ce);             
h23= a*b*(a-b)/(3*E*Ie*Ce);
h33= (a^3+b^3)/(3*E*Ie*Ce^2);
h11= h22;
h44= h33;
h14= -h23;
h41= h14;
h32= h23;

%Coeficientes de influencia, mancais flexivel e eixo rigido
h11_1= bet^2/(Kmy) + alp^2/(Kmy);
h14_1= bet/(Kmy*Ce) - alp/(Kmy*Ce);
h44_1= 1/(Kmy*Ce^2) + 1/(Kmy*Ce^2);

h22_1= bet^2/(Kmz) + alp^2/(Kmz);
h23_1= -bet/(Kmz*Ce) + alp/(Kmz*Ce);
h33_1= 1/(Kmz*Ce^2) + 1/(Kmz*Ce^2);

%Coeficientes de influencia, mancais flexivel e eixo flexivel
h11f= h11 + h11_1;
h14f= h14 + h14_1;
h44f= h44 + h44_1;
h22f= h22 + h22_1;
h23f= h23 + h23_1;
h33f= h33 + h33_1;

%Coeficientes de rigidez
delt1= h11f*h44f - h14f^2;
delt2= h22f*h33f - h23f^2;

k11= h44f/delt1;
k14= -h14f/delt1;
k44= h11f/delt1;
k22= h33f/delt2;
k23= -h23f/delt2;
k33= h22f/delt2;

K=[k11 0 0 k14;...
   0 k22 k23 0;...
   0 k23 k33 0;...
   k14 0 0 k44];

frq12=[];
frq22=[];
frq32=[];
frq42=[];
for i= 1:length(omg)
    
    GM= [G*omg(i) M;...
          M Z];
    
    KM= [K Z;...
         Z -M];
        
    res= eig(KM,-GM);
    omgn= sort(imag(res))/2/pi;
    
    frq12= [frq12 omgn(5)];
    frq22= [frq22 omgn(6)];
    frq32= [frq32 omgn(7)];
    frq42= [frq42 omgn(8)];
end

%Graficos caso 2
figure(4)
plot(omghz, omghz,'b--',omghz, omghz*2,'r--','linewidth', 1.5); hold on
plot(omghz, frq12,'m',omghz, frq22,'k',omghz, frq32,'g',omghz, frq42,'b','linewidth', 1.5)
grid on; grid minor;
set(gca,'fontsize',11)
legend({'\omega=\Omega','\omega=2\Omega','Mode 1 (BW)','Mode 1 (FW)','Mode 2 (BW)','Mode 2 (FW)'},'location','northwest')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Frequency \omega [Hz]')
title('Campbell diagram - off-centre disk, flexible shaft and bearings')

figure(5)
plot(omghz, omghz,'b--',omghz, omghz*2,'r--','linewidth', 1.5); hold on
plot(omghz, frq12,'m',omghz, frq22,'k',omghz, frq32,'g',omghz, frq42,'b','linewidth', 1.5)
grid on; grid minor; 
ylim([0 50]);
xlim([0 100]);
set(gca,'fontsize',11)
legend({'\omega=\Omega','\omega=2\Omega','Mode 1 (BW)','Mode 1 (FW)'},'location','northeast')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Frequency \omega [Hz]')
title('Campbell diagram - off-centre disk, flexible shaft and bearings')

figure(6)
plot(omghz, omghz,'b--',omghz, omghz*2,'r--','linewidth', 1.5); hold on
plot(omghz, frq12,'m',omghz, frq22,'k',omghz, frq32,'g',omghz, frq42,'b','linewidth', 1.5)
grid on; grid minor; 
ylim([0 500]);
xlim([0 200]);
set(gca,'fontsize',11)
legend({'\omega=\Omega','\omega=2\Omega','Mode 1 (BW)','Mode 1 (FW)'},'location','northwest')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Frequency \omega [Hz]')
title('Campbell diagram - off-centre disk, flexible shaft and bearings')

%%
%% Diferença entre resultados
figure(7)
plot(omghz, (frq12-frq1),'m',omghz, (frq22-frq2),'k',omghz, (frq32-frq3),'g',omghz, (frq42-frq4),'b','linewidth', 1.5)
grid on; grid minor;
ylim([-40 2])
set(gca,'fontsize',11)
legend({'Mode 1, F-R (BW)','Mode 1, F-R (FW)','Mode 2, F-R (BW)','Mode 2, F-R (FW)'},'location','southeast')
xlabel('Shaft speed \Omega [Hz]')
ylabel('Frequency difference \omega [Hz]')
title('Difference: flexible (F) minus rigid (R) bearings')



