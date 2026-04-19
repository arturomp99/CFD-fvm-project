function [A, b] = fvm_1D_euler( ...
        state, ...
        cells, ...
        t ...
    )
    %FVM_1D_EULER Discretización espacial de las ecuaciones de Euler 1D por volúmenes finitos
    %   Esta función coordina:
    %   1. Cálculo de flujos convectivos (via CONVECTIVE_FLUX_INTERPOLATOR)
    %   2. Aplicación de términos fuente (via SOURCE_TERMS)
    %   3. Ensamblado del sistema dw/dt = A*w + b
    %
    %   MODULARIDAD:
    %   ============
    %   El interpolador específico se selecciona en Config.CONVECTIVE_FLUX_INTERPOLATOR:
    %   - rusanov_interpolator
    %   - hllc_interpolator
    %   - linear_interpolator
    %
    %   Input
    %   ---------------------
    %   state : double (3*N×1)
    %       Vector de estado concatenado [ρ₁...ρₙ; (ρu)₁...(ρu)ₙ; E₁...Eₙ]
    %
    %   cells : struct array (1×N)
    %       Estructura de células procesadas por mesh_processor con:
    %       - .centroid: coordenadas de centroides
    %       - .connectivity: vecinos que comparten caras
    %   t : double
    %       Current time [s]    %
    %   Output:
    %   --------
    %   A : double/sparse (3*N×3*N)
    %       Matriz del operador espacial (típicamente =0 para esquemas no lineales)
    %
    %   b : double (3*N×1)
    %       Vector de términos independientes (divergencia de flujos + fuentes)

    [A, b] = convective_flux(state, cells);

    % Aplicar términos fuente si están habilitados
    centroids_x = get_cell_centroid_x(cells);
    [A_source, b_source] = Config.SOURCE_TERMS(state, centroids_x, t);

    A = A + A_source;
    b = b + b_source;
end
