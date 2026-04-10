function [w_new] = bw_euler( ...
        w, ...
        t, ...
        dt, ...
        f)
    %BW_EULER Integra un paso temporal la función f mediante un método de Euler
    %implícito de un paso y primer orden.
    %
    %   Inputs:
    %   -------
    %   w : array(double)
    %     Vector de estado en el instante t. Es un vector columna.
    %   t : double
    %     Instante temporal a partir del que integrar. [s]
    %   dt : double
    %     Paso temporal.  [s]
    %   f :
    %     Función que calcula la derivada temporal de w.
    %
    %   Outputs:
    %   --------
    %   w_new : Vector de estado en t + dt.

    [A, b] = f(w, t);
    I = eye(size(A));
    K = (I - A * dt);
    c = w + b * dt;
    w_new = linsolve(K, c);

end
