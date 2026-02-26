function [stop] = stop_at_floor(w, t)
%STOP_AT_FLOOR Documentación muy extensa e interesante.
%   Detalles varios del cálculo.
%
%   Inputs:
%   -------
%   w : Vector de estado.
%   t : Instante temporal. [s]
%
%   Outputs:
%   --------
%   stop : True si hay que parar de calcular.

    stop = true;
    if (w(2) >= 0.)
        stop = false;
    end

end
