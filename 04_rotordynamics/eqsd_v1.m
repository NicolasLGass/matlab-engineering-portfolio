function dx= eqsd_v1(t,x,omg,kd,cd,vars)

kyy= kd(1);
kyz= kd(2);
kzy= kd(3);
kzz= kd(4);

cyy= cd(1);
cyz= cd(2);
czy= cd(3);
czz= cd(4);

k= vars(1); 
c= vars(2); 
g= vars(3);
M= vars(4); 
me= vars(5); 

uy= x(1);
uz= x(2);
ym= x(3);
zm= x(4);
uyp= x(5);                  % uy ponto (derivada)
uzp= x(6);                  % uz ponto (derivada)

ymp= (-k*(zm-uz)-2*(kzy*ym+kzz*zm)+(czz/cyz)*(k*(ym-uy)+2*(kyy*ym+kyz*zm)))/ (2*czy-(2*czz*cyy)/cyz);
zmp= (-k*(ym-uy)-2*(kyy*ym+kyz*zm+cyy*ymp))/(2*cyz);

uypp= (me*(omg^2)*(cos(omg*t))-c*uyp-k*(uy-ym))/M;
uzpp= (me*(omg^2)*(sin(omg*t))-M*g-c*uzp-k*(uz-zm))/M;

dx= [uyp; uzp; ymp; zmp; uypp; uzpp];

end