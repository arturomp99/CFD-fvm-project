function [A, b] = convective_flux(state, cells, boundary_info)
    %CONVECTIVE_FLUX Despacha el cálculo del flujo convectivo al interpolador configurado.
    %   [A, b] = CONVECTIVE_FLUX(state, cells, boundary_info) delega el ensamblaje de las matrices de flujo convectivo
    %   a cualquier interpolador configurado en Config.CONVECTIVE_FLUX_INTERPOLATOR
    %
    %   Inputs:
    %   -------
    %   state : column vector (3*N x 1)
    %     Vector de estado concatenado [densidad; momento; energía total] para todas las N celdas.
    %   cells : struct array (1 x N)
    %     Contienen información de geometría y conectividad.
    %   boundary_info : struct (optional)
    %     Contiene la configuración de condiciones de contorno.
    %
    %   Outputs:
    %   --------
    %   A : matrix (3N x 3N)
    %     Matriz del operador de flujo convectivo.
    %   b : column vector (3N x 1)
    %     Independent terms vector.

    if nargin < 3
        boundary_info = struct('boundary_types', {{}});
    end

    [A, b] = Config.CONVECTIVE_FLUX_INTERPOLATOR(state, cells, boundary_info);
end
