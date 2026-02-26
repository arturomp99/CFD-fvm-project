function [A, b] = centered(V, S, c)
%UPWIND Calcula los flujos mediante el método centrado.
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

    for i=2:(n - 1)
        A(i, i) = c / V(i) / 2 * (S(i) - S(i + 1));
        A(i, i - 1) = c / V(i) / 2 * S(i);
        A(i, i + 1) = -c / V(i) / 2 * S(i + 1);
    end
    A(1, 1) = -c / V(1) / 2 * S(2);
    A(1, 2) = -c / V(1) / 2 * S(2);
    A(n, n - 1) = c / V(n) / 2 * S(n);
    A(n, n) = c / V(n) / 2 * S(n);
end
