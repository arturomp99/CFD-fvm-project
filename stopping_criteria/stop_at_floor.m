function [stop] = stop_at_floor(w, t)
    %STOP_AT_FLOOR Criterio de parada por tolerancia o convergencia
    %
    %   NOTA: Originalmente diseñada para un problema de proyectil, no para
    %   ecuaciones de Euler compresibles.
    %
    %
    %   Input
    %   ---------------------
    %   w : double (3*N×1)
    %       Vector de estado CFD [ρ; ρu; E] (mal interpretado actualmente)
    %
    %   t : double
    %       Tiempo actual [s] (no usado)
    %
    %   Output
    %   --------
    %   stop : logical
    %       Siempre true (detiene simulación inmediatamente)
    %

    stop = true;

    if (w(2) >= 0.)
        stop = false;
    end

end
