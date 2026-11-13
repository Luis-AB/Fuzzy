function [Xpl, Xpt, Div_max] = EEMeasure(X, dims, pop, Div_max)
%La forma de implementar es al final de cada iteración con
%las siguientes variable
% X = población
% dims = dimnesiones
% pop = tamaño de la población
% Div_max = diversidad maxima, la cual es un valor que se va actualizando
%           en el proceso iterativo y en la primer iteración es necesario
%           iniciarla en 0 en el script principal dentro de los parametros
%           iniciales como pop_size, dims, max_iter, etc.
% Xpl = porcentaje de exploración
% Xpt = porcentaje de explotación.


    Divj = zeros(1,dims);

    for d=1:dims
        mediana = median(X(:,d));
        for p=1:pop
            aux(p) = mediana - X(p,d);
        end
        Divj(d) = sum(aux)/pop;
    end

    Div = sum(Divj)/dims;

    if Div_max == 0
        Div_max = Div;
        Xpl = Div / Div_max * 100;

        Xpt = abs(Div - Div_max)/Div_max * 100;
    elseif Div > Div_max
        Div_max = Div;
        Xpl = Div / Div_max * 100;
        Xpt = abs(Div - Div_max)/Div_max * 100;
    else
        Xpl = Div / Div_max * 100;
        Xpt = abs(Div - Div_max)/Div_max * 100;
    end
end