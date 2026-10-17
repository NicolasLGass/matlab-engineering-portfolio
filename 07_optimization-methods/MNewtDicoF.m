%% Método do gradiente com busca dicotômica
% f- função avaliada
% x0- ponto inicial
% x1- ponto inicial (Dicotômica)
% x2- ponto final (Dicotômica)
% nim- número max de iterações
% errp- erro mínimo aceito
% itmaxDico- Iterações máximas do método dicotômico
% results- mostrar resultados dos métodos (0 ou 1)

function [val] = MNewtDicoF(f,x0,x1,x2,nim,errp,itmaxDico,results)
warning off;

val= x0;                                %Valor inicial
vecval= [x0];                             %Vetor de valores
[~, gi, Hi]= f([val]);
di= -inv(Hi)*gi;
% di= -inv(Hi)\gi;
alfa= MDicoF(f,x0,x1,x2,di,errp/10,errp,itmaxDico,0); %Dicotómica

for k=1:1:nim
    
    
    valalp= val+alfa*di;

% 	  falpp= f([valalp]);                 %Função em xi+alfa
%     fvalp= f([val]);                    %Função em xi
        
%     fprintf('k(for): %d\n', k);
%     fprintf('falp: %f | fval: %f \n', falpp, fvalp);
    val= valalp;

    vecval= [vecval, val];              %Atualiza vetor de valores
    
    [~, gt]= f([val]);
	dx= norm(gt,2);
    if errp > dx
        fprintf('Convergência satisfeita... \n');
        fprintf('Iterações: %d \n',k);
        break;
    end
    
    [~, gi, Hi]= f([val]);
    di= -inv(Hi)*gi;
    alfa= MDicoF(f,val,x1,x2,di,errp/10,errp,itmaxDico,0); %Dicotómica (a partir do ponto atual)
    
end %for k

if results==1
    
    figure(1)
    plot(vecval(1,:),vecval(2,:), 'bo'); hold on; grid on;
    plot(vecval(1,:),vecval(2,:), 'g');
    title('Iterates until convergence'); xlabel('x_1'); ylabel('x_2');

    if f(x0)<10
        faux= @(x1,x2) 10*(x2-x1^2)^2+(1-x1)^2;     %Função auxiliar plot
        zmax= 10;
    else
        faux= @(x1,x2) 200*(x2-x1^2)^2+(1-x1)^2;    %Função auxiliar plot
        zmax= 50;
    end
    
    for i=1:length(vecval)
        vecx(i)= vecval(1,i);
        vecy(i)= vecval(2,i);
        vecz(i)= f([vecval(1,i),vecval(2,i)]);
    end
    
    figure(4)
    [XXs, YYs] = meshgrid(linspace(-5, 5, 121));        % fsurf replacement (MATLAB/Octave)
    surf(XXs, YYs, arrayfun(faux, XXs, YYs), 'FaceAlpha', 0.5, 'EdgeColor', 'none'); hold on; grid on;
    plot3(vecval(1,1:end),vecval(2,1:end), vecz,'r', 'linewidth', 5);
    zlim([-2 zmax]); ylim([-1 2]); xlim([-2 2]);
    view(40,35);
end %if results
val;
end