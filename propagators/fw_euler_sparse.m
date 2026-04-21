function [w_new] = fw_euler_sparse(w, t, dt, f)
    %FW_EULER Integrador temporal explícito de Euler hacia adelante
    %
    %   Implementa el método de Euler explícito de primer orden para integrar
    %   el sistema de EDO: dw/dt = A*w + b
    %
    %   Input
    %   ---------------------
    %   w : double (3*N×1)
    %       Vector de estado en tiempo t [ρ; ρu; E]
    %
    %   t : double
    %       Tiempo actual [s] (pasado a la función problema)
    %
    %   dt : double
    %       Paso temporal [s] (debe satisfacer condición CFL)
    %
    %   f : function_handle
    %       Función problema que retorna [A, b] = f(w, t)
    %       donde dw/dt = A*w + b
    %
    %   Output
    %   --------
    %   w_new : double (3*N×1)
    %       Vector de estado en tiempo t + dt
    %
    %   Ver también: bw_euler, constant_dt, euler_cfl_dt

    [A, b] = f(w, t);
    I = speye(size(A));
    w_new = (I + A * dt) * w + b * dt;

end
