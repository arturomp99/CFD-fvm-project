function [w_new] = bw_euler_sparse( ...
        w, ...
        t, ...
        dt, ...
        f)
    %BW_EULER Integrador temporal implícito de Euler hacia atrás
    %
    %   Implementa el método de Euler implícito de primer orden para integrar
    %   el sistema de EDO: dw/dt = A*w + b
    %
    %   Input
    %   ---------------------
    %   w : double (3*N×1)
    %       Vector de estado en tiempo t [ρ; ρu; E]
    %
    %   t : double
    %       Tiempo actual [s]
    %
    %   dt : double
    %       Paso temporal [s] (sin restricción CFL)
    %
    %   f : function_handle
    %       Función problema que retorna [A, b] = f(w, t)
    %
    %   Output:
    %   --------
    %   w_new : double (3*N×1)
    %       Vector de estado en tiempo t + dt
    %

    [A, b] = f(w, t);

    I = speye(size(A));
    K = (I - A * dt);
    c = w + b * dt;
    w_new = K \ c;

end
