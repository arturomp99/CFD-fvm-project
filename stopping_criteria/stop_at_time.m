function [stop] = stop_at_time(t, t_end)
    %STOP_AT_TIME Criterio de parada basado en tiempo máximo de simulación
    %   detener cuando se alcanza un tiempo final especificado.
    %
    %   Input
    %   ---------------------
    %   t : double
    %       Tiempo actual de simulación [s]
    %
    %   t_end : double
    %       Tiempo final deseado [s]
    %
    %   Output
    %   --------
    %   stop : logical
    %       true si t >= t_end, false en caso contrario

    if (t >= t_end)
        stop = true;
    else
        stop = false;
    end

end
