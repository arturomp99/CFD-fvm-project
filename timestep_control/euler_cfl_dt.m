function dt = euler_cfl_dt(w, t, cells, CFL)
    %EULER_CFL_DT Calculadora de paso temporal adaptativo basado en condición CFL
    %
    %   Calcula el paso temporal máximo permitido por la condición de estabilidad
    %   CFL (Courant-Friedrichs-Lewy) para las ecuaciones de Euler compresibles.
    %   La condición CFL asegura que la información no se propague
    %   más rápido de lo que el esquema numérico puede resolverla.
    %
    %   Input
    %   ---------------------
    %   w : double (3*N×1)
    %       Vector de estado [ρ; ρu; E]
    %       
    %   t : double
    %       Tiempo actual [s] (no usado, para compatibilidad)
    %       
    %   cells : struct array (1×N)
    %       Estructura de células con campo .centroid
    %       
    %   CFL : double
    %       Número de Courant deseado (0 < CFL ≤ 1)
    %
    %   Output
    %   --------
    %   dt : double
    %       Paso temporal máximo permitido [s]
    %       

    num_cells = length(cells);

    % Extracción de variables conservadas del vector de estado
    rho = w(1:num_cells);                           % Densidad [kg/m³]
    rhou = w(num_cells + 1:2 * num_cells);          % Momentum [kg/(m²·s)]
    E = w(2 * num_cells + 1:3 * num_cells);         % Energía total [J/m³]

    % Limitación para evitar valores no físicos (densidad ~0)
    rho = max(rho, 1e-10);

    % Cálculo de variables primitivas
    u = rhou ./ rho;                                % Velocidad [m/s]
    p = (Air.GAMMA - 1) .* (E - 0.5 .* (rhou .^ 2) ./ rho); % Presión [Pa]
    
    % Limitación para evitar presiones negativas
    p = max(p, 1e-10);

    % Velocidad del sonido
    c = sqrt(Air.GAMMA .* p ./ rho);                % [m/s]
    
    % Velocidades características máximas
    lambda = abs(u) + c;                            % |u| + c [m/s]

    % Extracción de coordenadas x de centroides de células
    x = get_cell_centroid_x(cells);

    % CRÍTICO: eliminar coordenadas x repetidas 
    % Esto puede ocurrir en mallas 2D donde varias células tienen la misma x
    x_unique = unique(sort(x));

    if numel(x_unique) < 2
        error('Not enough distinct x-coordinates to compute CFL timestep.');
    end

    % Cálculo del tamaño mínimo de celda Δx_min
    dx_min = min(diff(x_unique));
    dx_min = max(dx_min, 1e-6);  % Evitar Δx demasiado pequeño

    % Velocidad característica máxima en todo el dominio
    max_lambda = max(lambda);
    max_lambda = max(max_lambda, 1e-10);  % Evitar división por cero

    % Aplicación de la condición CFL: dt = CFL * Δx_min / λ_max
    dt = CFL * dx_min / max_lambda;
end
