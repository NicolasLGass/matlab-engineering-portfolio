%Função f4
function [f, g, H] = func4(x)

A4= [1 1 1 1; ...
     1 3 3 3; ...
     1 3 5 5; ...
     1 3 5 7];
b4= [10; 28; 42; 50];
f= (0.5*x'*A4*x)-(b4'*x);

g= A4*x-b4;         %Gradiente da função

H= A4;              %Hessiana da função
end