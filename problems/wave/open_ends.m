function [A, b] = open_ends(V, S, c)
%OPEN_ENDS Calcula flujos debido a extremos abiertos.
%
%   Inputs:
%   -------
%   V :
%     Vector de volúmenes de cada celda.
%   S :
%     Vector de áreas entre celdas.
%   c :
%     Velocidad de propagación. [m / s]
%
%   Outputs:
%   --------
%   A :
%     Matriz jacobiana de dw_dt evaluada en el instante t.
%   b :
%     Vector de términos independientes evaluado en el instante t.
    
    n = length(V);
    A = zeros(n, n);
    b = zeros(n, 1);

    if (c > 0)
        % Sale flujo de la celda n.
            A(n, n) = - c * S(n + 1) / V(n);
    else
        % Sale flujo de la celda 1.
        A(1, 1) = - c * S(1) / V(1);
    end
end
