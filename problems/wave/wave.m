function [A, b] = wave(w, t, V, S, c, fluxes_method, bcs)
%WAVE Calcula la derivada temporal del vector de estado en una onda que
%viaja en una sola dirección.
%
%   Inputs:
%   -------
%   w :
%     Valor medio del vector de estado en cada celda. Es un vector columna.
%   t :
%     Tiempo actual. [s]
%   V :
%     Vector de volúmenes de cada celda.
%   S :
%     Vector de áreas entre celdas.
%   c :
%     Velocidad de propagación. [m / s]
%   fluxes_method : Method used to compute the fluxes.
%   bcs : Vector of boundary conditions.
%
%   Outputs:
%   --------
%   A :
%     Matriz jacobiana de dw_dt evaluada en el instante t.
%   b :
%     Vector de términos independientes evaluado en el instante t. Es un
%     vector columna.

    [A1, b1] = fluxes_method(V, S, c);
    [A2, b2] = bcs(V, S, c, t);
    A = A1 + A2;
    b = b1 + b2;

end

