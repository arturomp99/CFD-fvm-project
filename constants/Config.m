classdef Config
    %CONFIG Configuración central del solver CFD de volúmenes finitos
    %
    %   Esta clase contiene todos los parámetros configurables del solver
    %   para personalizar simulaciones sin modificar el código principal.
    %   
    %   🕒 CONTROL TEMPORAL:
    %   ===================

    properties (Constant)
        % Tiempo inicial y final de simulación [s]
        T0 = 0.0;         % Tiempo inicial (típicamente 0)
        T_END = 0.2;      % Tiempo final (ajustar según el problema)
        
        % Intervalo de muestreo para resultados [s]
        SAMPLE_DT = 1e-3; % Frecuencia de guardado (menor = más puntos temporales)

        % 🏗️ CONFIGURACIÓN DE CONDICIONES DE FRONTERA:
        % ============================================
        % 
        % Tipos disponibles por superficie:
        % 'open'     - Frontera transmisiva (ondas pasan libremente)
        % 'wall'     - Pared sólida reflexiva (velocidad normal = 0)
        % 'velocity' - Velocidad impuesta (requiere BOUNDARY_VELOCITIES)
        %
        % ⚠️ ORDEN IMPORTANTE: Debe coincidir con FilePaths.BOUNDARY_CONDITIONS
        BOUNDARY_TYPES = { ...
                              'open', ...     % 1: bc_bottom (entrada/salida libre)
                              'velocity', ... % 2: bc_left (velocidad prescrita)
                              'open', ...     % 3: bc_right (salida libre)
                              'open', ...     % 4: bc_top (frontera abierta)
                              'open' ...      % 5: bc_whole_contour (contorno completo)
                          };

        % Valores de velocidad para condiciones 'velocity' [m/s]
        % Solo se usan cuando el BOUNDARY_TYPES correspondiente es 'velocity'
        % Usar NaN para superficies que no requieren velocidad específica
        BOUNDARY_VELOCITIES = [NaN, 100.0, NaN, NaN, NaN];

        % 🎯 DEFINICIÓN DEL PROBLEMA:
        % ===========================
        
        % Condiciones iniciales - función que define el estado inicial
        % Opciones predefinidas: uniform(), sod()
        % Formato: @(pos) función_condicion_inicial(parámetros..., pos)
        INITIAL_CONDITIONS = @(pos) ... 
            uniform( ...
                Air.SEA_LEVEL_PRESSURE, ... % 101325 Pa (presión atmosférica)
                Air.SEA_LEVEL_DENSITY, ...  % 1.225 kg/m³ (densidad del aire)
                0.0, ...                    % 0 m/s (velocidad inicial)
                pos ...                     % coordenadas x de centroides
            );
            
        % Ejemplo alternativo - Tubo de choque de Sod (descomentado para usar):
        %     sod( ...
        %         struct('left', 1.0, 'right', 0.1), ...    % presión [Pa]
        %         struct('left', 1.0, 'right', 0.125), ...  % densidad [kg/m³]
        %         struct('left', 0.0, 'right', 0.0), ...    % velocidad [m/s]
        %         0.5, ...                                  % posición discontinuidad [m]
        %         pos ...
        %     );

        % Términos fuente - función para fuentes/sumideros de masa, momentum o energía
        % Actualmente retorna cero (sin fuentes). Modificar para casos específicos.
        SOURCE_TERMS = @(state, pos, t) ...
            point_source(); % Función stub - retorna [0; 0; 0]

        % Condiciones de frontera adicionales - para BCs complejas dependientes del tiempo
        % Actualmente no usado (retorna 0). Extender para casos avanzados.
        BOUNDARY_CONDITIONS = @(state, pos, t) ...
            0;

        % 🔧 CONFIGURACIÓN DEL SOLVER:
        % ============================

        % Control de implicitidad para términos fuente
        IS_SOURCE_IMPLICIT = false; % true = implícito (mayor estabilidad), false = explícito

        % Control de implicitidad para términos fuente
        IS_SOURCE_IMPLICIT = false; % true = implícito (mayor estabilidad), false = explícito

        % Esquema de interpolación para flujos convectivos
        % Opciones disponibles en convective_flux/interpolators/:
        % - hllc_interpolator: HLLC Riemann solver (alta resolución, captura choques)
        % - rusanov_interpolator: 1er orden, estable, difusivo (más robusto)
        % - upwind_interpolator: 1er orden, estable (similar a Rusanov)  
        % - linear_interpolator: 2do orden, menos difusivo, puede ser inestable
        CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) ... 
            hllc_interpolator(state, cells); % Esquema HLLC por defecto

        % Integrador temporal - define método de avance en tiempo
        % - fw_euler: Explícito (rápido, limitado por CFL)
        % - bw_euler: Implícito (estable, requiere solver lineal)
        PROPAGATOR = @(state, time, d_time, problem) ...
            bw_euler(state, time, d_time, problem); % Implícito por defecto

        % Calculadora de paso temporal
        % - constant_dt: Paso fijo (simple, puede ser ineficiente)
        % - euler_cfl_dt: Adaptativo basado en condición CFL (recomendado)
        TIMESTEP_CALCULATOR = @(w, t) ...
            constant_dt(w, t, 5e-5); % dt = 5×10⁻⁵ s (ajustar según estabilidad)

        % Condición de parada de la simulación
        % - stop_at_time: Para en tiempo específico
        % - stop_at_floor: Para cuando cambios < tolerancia
        STOPPING_CONDITION = @(w, t) ...
            stop_at_time(t, Config.T_END);

        % Gestor de resultados - controla cómo se almacenan/muestrean los datos
        % - sample_results: Muestrea cada SAMPLE_DT segundos (recomendado)
        % - save_last_results: Solo guarda el resultado final
        % - discard_results: No guarda resultados (solo para análisis de rendimiento)
        RESULTS_MANAGER = @(w, t, old_results) ...
            sample_results(w, t, old_results, Config.SAMPLE_DT);

        % 📊 CONFIGURACIÓN DE VISUALIZACIÓN:
        % =================================
        
        % Control de suavizado en gráficos espacio-tiempo
        % false = imagesc (píxeles exactos, más rápido)
        % true = contourf (contornos suavizados, más estético)
        DIFFUMINATED_VISUALIZATION = false;
        
    end
    %
    % 📖 GUÍA DE CONFIGURACIÓN RÁPIDA:
    % ================================
    % 
    % Para un CASO ESTÁNDAR (flujo uniforme):
    % - Mantener configuración actual
    % 
    % Para el TUBO DE CHOQUE DE SOD:
    % - Descomentar la sección sod() en INITIAL_CONDITIONS
    % - Usar BOUNDARY_TYPES = {'open', 'open', ...} en ambos extremos
    % 
    % Para MAYOR PRECISIÓN:
    % - Mantener hllc_interpolator (ya configurado)
    % - Reducir TIMESTEP_CALCULATOR a dt más pequeño
    % - Usar DIFFUMINATED_VISUALIZATION = true
    % 
    % Para MAYOR ESTABILIDAD:
    % - Cambiar a rusanov_interpolator (más difusivo pero estable)
    % - Mantener bw_euler (implícito)
    % - Usar euler_cfl_dt para paso adaptativo
    %
    % Ver README.md para guías detalladas y ejemplos de uso.
    %

end
