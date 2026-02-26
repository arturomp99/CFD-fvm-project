function [A, b] = cannon(w, t, a_x, a_y)
%CANNON Calcula la derivada temporal del vector de estado en un tiro
%parabólico.
%
%   Inputs:
%   -------
%   w : 
%     Vector de estado, formado por la posición X, la posición Y, la
%     velocidad según X y la velocidad según Y. Es un vector columna.
%     [[m], [m], [m/s], [m/s]]
%   t : double
%     Tiempo actual. [s]
%   a_x : double
%     Aceleración según x. [m / s ** 2]
%   a_y : Aceleración según y: [m / s ** 2]
%
%   Outputs:
%   --------
%   A : 
%     Matriz jacobiana de dw_dt evaluada en el instante t.
%   b :
%     Vector de términos independientes evaluado en el instante t. Es un
%     vector columna

%    dw_dt = zeros(size(w));
%    dw_dt(1) = w(3);
%    dw_dt(2) = w(4);
%    dw_dt(3) = a_x;
%    dw_dt(4) = a_y;
%    A : matriz jacobiana de dw_dt;
    
    A = [[0, 0, 1, 0];
        [0, 0, 0, 1]; 
        [0, 0, 0, 0]; 
        [0, 0, 0, 0]];
    b = [0, 0, a_x, a_y]';

end
