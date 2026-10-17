%Função f1
function [f, g, H] = func2(x)

f= 200*(x(2)-x(1)^2)^2+(1-x(1))^2;
g= [(2*x(1) - 800*x(1)*(- x(1)^2 + x(2)) - 2); (- 400*x(1)^2 + 400*x(2))];       %Gradiente da função
H= [ 2400*x(1)^2 - 800*x(2) + 2, -800*x(1); -800*x(1),    400];                  %Hessiana da função

end