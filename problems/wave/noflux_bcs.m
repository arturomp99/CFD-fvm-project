function [A, b] = noflux_bcs(n)
%NOFLUX_BCS Calcula flujos nulos en las condiciones de contorno.
%
%   Inputs:
%   -------
%   n : Número de celdas.
%
%   Outputs:
%   --------
%   A :
%     Matriz jacobiana de dw_dt evaluada en el instante t.
%   b :
%     Vector de términos independientes evaluado en el instante t.
    
    A = zeros(n, n);
    b = zeros(n, 1);

end
