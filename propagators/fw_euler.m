function [w_new] = fw_euler(w, t, dt, f)
%FW_EULER Integra un paso temporal la función f mediante un método de Euler
%explícito de un paso y primer orden.
%
%   Inputs:
%   -------
%   w :
%     Vector de estado en el instante t. Es un vector columna.
%   t :
%     Instante temporal a partir del que integrar. [s]
%   dt :
%     Paso temporal. [s]
%   f :
%     Función que calcula la derivada temporal de w.
%
%   Outputs:
%   --------
%   w_new :
%     Vector de estado en t + dt. Es un vector columna.

    [A, b] = f(w, t);
    I = eye(size(A));
    w_new = (I + A * dt) * w + b * dt;

end

