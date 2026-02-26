function [dt] = CFL_dt(CFL, c, dx)
%CFL_DT Documentación muy extensa e interesante.
%   Detalles varios del cálculo.
%
%   Inputs:
%   -------
%   w :
%     Vector de estado.
%   t :
%     Instante temporal. [s]
%   CFL :
%     Courant-Friedrichs-Levvy constant.
%   c :
%     Propagation speed. [m / s]
%   dx :
%     Mesh size vector. [m]
%
%   Outputs:
%   --------
%   dt :
%     Paso temporal a dar. [s]

    dt = CFL * min(dx ./ c);

end
