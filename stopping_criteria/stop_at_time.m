function [stop] = stop_at_time(t, t_end)
    %STOP_AT_TIME Documentación muy extensa e interesante.
    %   Detalles varios del cálculo.
    %
    %   Inputs:
    %   -------
    %   t : Instante temporal. [s]
    %   t_end : Tiempo en el que parar. [s]
    %
    %   Outputs:
    %   --------
    %   stop : True si hay que parar de calcular.

    if (t >= t_end)
        stop = true;
    else
        stop = false;
    end

end
