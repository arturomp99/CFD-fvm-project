classdef Config
    %CONFIG Configuración central del solver CFD de volúmenes finitos
    %

    properties (Constant)
        T0 = 0.0; % Tiempo inicial (típicamente 0)
        T_END = 0.2; % Tiempo final (ajustar según el problema)
        SAMPLE_DT = 1e-3; % Intervalo de muestreo para resultados [s]

        % ============================================
        % CONFIGURACIÓN DE CONDICIONES DE FRONTERA:
        % ============================================
        %
        % Tipos disponibles por superficie:
        %   'open'     - Frontera transmisiva (ondas pasan libremente)
        %   'wall'     - Pared sólida reflexiva (velocidad normal = 0)
        %   'velocity' - Velocidad impuesta (requiere BOUNDARY_VELOCITIES)
        % ORDEN IMPORTANTE: Debe coincidir con FilePaths.BOUNDARY_CONDITIONS
        BOUNDARY_TYPES = { ...
                              'open', ... % 1: bc_bottom (entrada/salida libre)
                              'velocity', ... % 2: bc_left (velocidad prescrita)
                              'open', ... % 3: bc_right (salida libre)
                              'open', ... % 4: bc_top (frontera abierta)
                              'open' ... % 5: bc_whole_contour (contorno completo)
                          };

        % Valores de velocidad para condiciones 'velocity' [m/s]
        % Solo se usan cuando el BOUNDARY_TYPES correspondiente es 'velocity'
        % Usar NaN para superficies que no requieren velocidad específica
        BOUNDARY_VELOCITIES = [ ...
                                   NaN, ...
                                   100.0, ...
                                   NaN, ...
                                   NaN, ...
                                   NaN ...
                               ];

        % ============================
        % DEFINICIÓN DEL PROBLEMA:
        % ===========================
        % Condiciones iniciales - función que define el estado inicial
        % Opciones predefinidas: uniform(), sod()
        INITIAL_CONDITIONS = @(pos) ... % funciónn de la posición
            uniform( ...
            Air.SEA_LEVEL_PRESSURE, ... % 101325 Pa (presión atmosférica)
            Air.SEA_LEVEL_DENSITY, ... % 1.225 kg/m³ (densidad del aire)
            0.0, ... % 0 m/s (velocidad inicial)
            pos ... % coordenadas x de centroides
        );

        % Ejemplo alternativo - Tubo de choque de Sod (descomentado para usar):
        %     sod( ...
        %         struct('left', 1.0, 'right', 0.1), ...    % presión [Pa]
        %         struct('left', 1.0, 'right', 0.125), ...  % densidad [kg/m³]
        %         struct('left', 0.0, 'right', 0.0), ...    % velocidad [m/s]
        %         0.5, ...                                  % posición discontinuidad [m]
        %         pos ...
        %     );

        % Términos fuente - función para fuentes/sumideros de masa, momento o energía
        SOURCE_TERMS = @(state, pos, t) ...
            point_source();

        % Condiciones de frontera adicionales - para BCs complejas dependientes del tiempo
        % Actualmente no usado (retorna 0). Extender para casos avanzados.
        BOUNDARY_CONDITIONS = @(state, pos, t) ...
            0;

        % ============================
        % CONFIGURACIÓN DEL SOLVER:
        % ============================
        IS_SOURCE_IMPLICIT = false; % Control de implicitidad para términos fuente

        % Esquema de interpolación para flujos convectivos
        % Opciones disponibles en convective_flux/interpolators/:
        %   - hllc_interpolator: HLLC Riemann solver (alta resolución, captura choques)
        %   - rusanov_interpolator: 1er orden, estable, difusivo (más robusto)
        %   - upwind_interpolator: 1er orden, estable (similar a Rusanov)
        %   - linear_interpolator: 2do orden, menos difusivo, puede ser inestable
        CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) ...
            hllc_interpolator(state, cells); % Esquema HLLC por defecto

        % Integrador temporal - define método de avance en tiempo
        %   - fw_euler: Explícito (rápido, limitado por CFL)
        %   - bw_euler: Implícito (estable, requiere solver lineal)
        PROPAGATOR = @(state, time, d_time, problem) ...
            bw_euler(state, time, d_time, problem); % Implícito por defecto

        % Calculadora de paso temporal
        %   - constant_dt: Paso fijo (simple, puede ser ineficiente)
        %   - euler_cfl_dt: Adaptativo basado en condición CFL (recomendado)
        TIMESTEP_CALCULATOR = @(w, t) ...
            constant_dt(w, t, 5e-5); % dt = 5×10⁻⁵ s (ajustar según estabilidad)

        % Condición de parada de la simulación
        %   - stop_at_time: Para en tiempo específico
        %   - stop_at_floor: Para cuando cambios < tolerancia
        STOPPING_CONDITION = @(w, t) ...
            stop_at_time(t, Config.T_END);

        % Gestor de resultados
        %   - sample_results: Muestrea cada SAMPLE_DT segundos (recomendado)
        %   - save_last_results: Solo guarda el resultado final
        %   - discard_results: No guarda resultados (solo para análisis de rendimiento)
        RESULTS_MANAGER = @(w, t, old_results) ...
            sample_results(w, t, old_results, Config.SAMPLE_DT);

        % =================================
        % CONFIGURACIÓN DE VISUALIZACIÓN:
        % =================================

        % Control de suavizado en gráficos espacio-tiempo
        DIFFUMINATED_VISUALIZATION = false;

    end

    % Ver README.md para guía de uso detalladas

end
