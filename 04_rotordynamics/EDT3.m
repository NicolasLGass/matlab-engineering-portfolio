function Eqs1 = EDT3(t,x)

T= 0.075;                   %Torque aplicado

%Parametros do material

E= 200*(10^9);              %Modulo de elasticidade [Pa]
rho= 7850;                  %Densidade [kg/m³]
Ps= 0.3;                    %Poisson
me= 10^(-3);                %Desbalanceamento [kg/m]

%Eixo
De= 0.01;                   %Diametro eixo [m]
Ce= 0.80;                   %Comprimento eixo [m]
Ae= pi*(De/2)^2;            %Area do eixo [m²]
% K= (Ae*E)/Ce;               %Rigidez do eixo
I= (pi*(De/2)^4)/4;         %Momento de inércia do eixo [m^4]
K= (48*E*I)/(Ce^3);         %Rigidez do eixo
c= 10^(-3)*K;               %Coeficiente de amortecimento [Ns/m]

%Disco
Dd= 0.09;                   %Diametro disco [m]
Cd= 0.02;                   %Comprimento disco [m]
Ad= pi*(Dd/2)^2;            %Area do disco [m²]
Vd= Ad*Cd;                  %Volume do disco [m³]
md= Vd*rho;                 %Massa do disco []
Jg= (md*(Dd^2))/8;          %Momento de inercia [m^4]

Ex= me/md;                  %Excentricidade

wn= sqrt(K/md)*(60/(2*pi)); %Freq natural [rpm]

% Equacoes do movimento e variaveis

uy= x(1);
uz= x(2);
theta= x(3);
uyp= x(4);
uzp= x(5);
thetap= x(6);

thetapp= (T+K*uz*Ex*cos(theta)-K*uy*(Ex*sin(theta)) - c*uyp*(Ex*sin(theta))+c*uzp*(Ex*cos(theta)))/Jg;
uypp= (me*(thetapp*sin(theta)+(thetap^2)*cos(theta))-c*uyp-K*uy)/md;
uzpp= (me*(-thetapp*cos(theta)+(thetap^2)*sin(theta))-md*9.81-c*uzp-K*uz)/md;

% ddtheta = (T+k*uz*e*cos(theta)-k*uy*e*sin(theta) - c*duy*e*sin(theta)+c*duz*e*cos(theta))/Jg;
% dduy = (me*(ddtheta*sin(theta)+(dtheta^2)*cos(theta))-c*duy-k*uy)/M;
% dduz = (me*(-ddtheta*cos(theta)+(dtheta^2)*sin(theta))-M*g-c*duz-k*uz)/M;
    

Eqs1= [uyp; uzp; thetap; uypp; uzpp; thetapp];

end

