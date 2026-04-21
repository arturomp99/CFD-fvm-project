% ===========================================================================
% MAIN.M - Solver CFD de Volúmenes Finitos para Ecuaciones de Euler 1D
% ===========================================================================
%
% Este script ejecuta una simulación completa CFD usando el método de
% volúmenes finitos para resolver las ecuaciones de Euler compresibles 1D.
%
% FLUJO DE EJECUCIÓN:
% ==================
% 1. Configuración inicial y paths
% 2. Procesamiento de malla y geometría
% 3. Configuración del problema físico
% 4. Resolución temporal del sistema
% 5. Post-procesamiento y visualización
%
% Ver Config.m para personalizar parámetros de simulación.

clc; clear all; close all;
% ===========================================================================
% 1. CONFIGURACIÓN DE PATHS Y CARGA DE MÓDULOS
% ===========================================================================
% Agrega todos los directorios necesarios al path de MATLAB para acceder
% a las funciones del solver CFD

addpath('mesh_processing'); % Procesamiento de malla y geometría
addpath(genpath('utils')); % Utilidades (termodinámica, fronteras, etc.)
addpath('constants'); % Constantes físicas y configuración
addpath('initial_conditions'); % Condiciones iniciales predefinidas
addpath('problems/fvm_1D_euler'); % Implementación específica de Euler 1D
addpath(genpath('convective_flux')); % Esquemas de interpolación de flujos
addpath('stopping_criteria'); % Criterios de parada de simulación
addpath('timestep_control'); % Control de paso temporal
addpath('propagators'); % Integradores temporales
addpath('results_manager'); % Gestión de resultados y muestreo
addpath('sources'); % Términos fuente (actualmente deshabilitados)
addpath('visualizer'); % Herramientas de visualización
addpath('console_logs'); % Mensajes informativos
addpath('boundary_condition'); % Condiciones de frontera

%% ========================================================================
%% 2. PROCESAMIENTO DE MALLA
%% ========================================================================
% Carga archivos de malla (.dat) y construye estructuras de datos optimizadas
% para cálculos FVM con información geométrica completa

welcome_msg(); % Mensaje de bienvenida con información del proyecto

% Rutas a archivos de malla (configuradas en FilePaths.m)
nodes_file = FilePaths.NODES; % Coordenadas de nodos
cells_file = FilePaths.CELLS; % Conectividad de células
bc_files = FilePaths.BOUNDARY_CONDITIONS; % Superficies de frontera

starting_mesh_processing_msg(); % Mensaje informativo con timer
cells = mesh_processor(nodes_file, cells_file, bc_files);
finishing_mesh_processing_msg(); % Tiempo transcurrido

% ESTRUCTURA DE DATOS RESULTANTE:
% cells(i) = estructura con geometría completa de la célula i
%   .nodes              - Coordenadas [x, y] de vértices de la célula [m]
%   .faces              - Definiciones de caras como endpoints [m]
%   .volume             - Área de la célula (en 2D) [m²]
%   .centroid           - Centro geométrico [cx, cy] [m]
%   .area               - Longitudes de aristas (áreas de caras en 2D) [m]
%   .normals            - Vectores normales unitarios salientes por cara
%   .connectivity       - Índices de células vecinas que comparten caras
%

%% ========================================================================
%% 3. CONFIGURACIÓN DEL PROBLEMA FÍSICO
%% ========================================================================
% Define condiciones iniciales, esquemas numéricos y parámetros de simulación

% VECTOR DE ESTADO CONCATENADO:
% El solver utiliza un vector de estado agrupado por variable de tamaño 3*N×1:
% w = [ρ₁; ρ₂; ...; ρₙ;           % Densidad para cada celda [kg/m³]
%      (ρu)₁; (ρu)₂; ...; (ρu)ₙ;  % Momentum para cada celda [kg/(m²·s)]
%      E₁; E₂; ...; Eₙ]           % Energía total para cada celda [J/m³]

% Configuración de condiciones iniciales desde Config.m
initial_conditions = Config.INITIAL_CONDITIONS;

% Extrae coordenadas x de centroides para condiciones iniciales espacialmente variables
centroids_x = get_cell_centroid_x(cells);

% Calcula vector de estado inicial aplicando condición inicial configurada
w0 = initial_conditions(centroids_x);

% FUNCIÓN DEL PROBLEMA:
% Define la función que calcula [A, b] tal que dw/dt = A*w + b
% Nota: Usando versión implícita que puede manejar matrices no-cero A
problem = @(state, time) fvm_1D_euler_implicit(state, cells, time);

% COMPONENTES DEL SOLVER CONFIGURABLES:
% Todos definidos en Config.m para fácil personalización

% Integrador temporal (explícito vs implícito)
propagator = Config.PROPAGATOR;

% Control de paso temporal (fijo vs adaptativo CFL)
timestep_calculator = Config.TIMESTEP_CALCULATOR;

% Criterio de parada (tiempo vs convergencia)
stopping_condition = Config.STOPPING_CONDITION;

% Gestor de resultados - muestreo temporal para post-procesamiento
% Configurado para muestrear cada SAMPLE_DT segundos y crear matriz
% estructurada compatible con visualizer()
manager = Config.RESULTS_MANAGER;

%% ========================================================================
%% 4. RESOLUCIÓN TEMPORAL
%% ========================================================================
% Ejecuta bucle de integración temporal hasta satisfacer criterio de parada

starting_solver_msg(); % Mensaje informativo con timer

results = solver( ...
    w0, ... % Vector de estado inicial [3*N×1]
    Config.T0, ... % Tiempo inicial [s]
    problem, ... % Función que calcula dw/dt = A*w + b
    propagator, ... % Integrador temporal (fw_euler/bw_euler)
    timestep_calculator, ... % Calculadora de dt (constante/CFL adaptativo)
    stopping_condition, ... % Criterio de parada (tiempo/convergencia)
    manager ... % Gestor de muestreo de resultados
);

finishing_solver_msg(); % Tiempo transcurrido total

% FORMATO DE RESULTADOS:
% results es una matriz (num_samples × (3*N + 1)) donde:
%   results(i, 1)         → tiempo t [s] en muestra i
%   results(i, 2:N+1)     → densidad ρ para cada celda [kg/m³]
%   results(i, N+2:2N+1)  → momentum ρu para cada celda [kg/(m²·s)]
%   results(i, 2N+2:3N+1) → energía total E para cada celda [J/m³]

%% ========================================================================
%% 5. POST-PROCESAMIENTO Y VISUALIZACIÓN
%% ========================================================================
% Genera gráficos espacio-tiempo de todas las variables termodinámicas

visualizer(results, centroids_x);
