function [xInt, yInt]= intrpElipse(a, b)
x0=0;   % Centro da elipse em 0
y0=0;

r= linspace(-pi,pi,60);
x= x0+a*cos(r);
y= y0+b*sin(r);

dist = sqrt(diff(x).^2+diff(y).^2);

%Reparametrizando
tm = cumsum([0 dist]);
r = linspace(0,tm(end));
xInt = interp1(tm,x,r,"spline");
yInt = interp1(tm,y,r,"spline");
end