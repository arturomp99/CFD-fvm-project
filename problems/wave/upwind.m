function [A, b] = upwind(V, S, c)
%UPWIND Calcula los flujos mediante el método upwind.
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
        for i=1:(n - 1)
            % Sale flujo de la celda i hacia la celda i + 1. El vector de
            % estado de la celda i pierde algo de valor.
            A(i, i) = - c * S(i + 1) / V(i);
            % Entra flujo en la celda i + 1 desde la celda i. El vector
            % de estado de la celda i + 1 gana algo de valor.
            A(i + 1, i) = c * S(i + 1) / V(i + 1);
        end
    else
        for i=2:n
            % Sale flujo de la celda i hacia la celda i - 1. El vector de
            % estado de la celda i pierde algo de valor.
            A(i,  i) = - c * S(i) / V(i);
            % Entra flujo en la celda i - 1 desde la celda i. El vector
            % de estado de la celda i - 1 gana algo de valor.
            A(i - 1, i) = c * S(i) / V(i - 1);
        end
    end
end
