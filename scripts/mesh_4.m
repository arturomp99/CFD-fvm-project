%% CREATE_MISSING_REPORT_FIGURES_V7
% Genera las figuras principales para el informe del proyecto de Euler 1D.
% Versión revisada:
% - Coste computacional frente al número de celdas
% - Referencia de orden p=1 para el caso de Sod
% - GCI y rango asintótico calculados sobre TODOS los tripletes consecutivos
% - Gráficas GCI/AR en formato línea, no barras
% - Figura adicional del orden del coste computacional
% - En verificación y comparación se usa energía interna específica
 
clc; clear; close all;
 
%% =========================
%  SETTINGS
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

% ==========================
SETTINGS.PROJECT_ROOT = fileparts(mfilename('fullpath'));
SETTINGS.SAVE_DIR = fullfile(SETTINGS.PROJECT_ROOT, 'figures_report_auto_v6');
SETTINGS.RUN_IMPLICIT = true;
SETTINGS.SAVE_FIGURES = true;
SETTINGS.SAVE_MAT = true;
 
SETTINGS.RUNTIME_REPEATS = 1;         % repetir tiempos para reducir ruido
SETTINGS.RUNTIME_DROP_FIRST = true;   % descarta primera corrida (warm-up)
SETTINGS.COST_ORDER_MIN_CELLS = 258;  % usar mallas medias-finas para el orden del coste
 
SETTINGS.MESH_LIST = [66, 130, 258, 514, 1026];
SETTINGS.DT_FINEST = 1e-5;
SETTINGS.T0 = 0.0;
SETTINGS.T_END = 0.2;
SETTINGS.X0 = 0.5;
SETTINGS.GAMMA = Air.GAMMA;
 
SETTINGS.VERIFICATION_MESH = 1026;
SETTINGS.COMPARISON_MESH = 514;
SETTINGS.REFINEMENT_PLOT_SKIP = [1026];
SETTINGS.ORDER_REFERENCE = 1.0;   % Para Sod con discontinuidades
 
SETTINGS.IMPLICIT_MESH_LIST = [66, 130, 258];
SETTINGS.IMPLICIT_DT = [2e-3, 1e-3, 5e-4];
 
SETTINGS.LEFT.rho = 1.0;
SETTINGS.LEFT.u   = 0.0;
SETTINGS.LEFT.p   = 1.0;
SETTINGS.RIGHT.rho = 0.125;
SETTINGS.RIGHT.u   = 0.0;
SETTINGS.RIGHT.p   = 0.1;
 
if ~exist(SETTINGS.SAVE_DIR, 'dir')
    mkdir(SETTINGS.SAVE_DIR);
end
 
setup_project_paths(SETTINGS.PROJECT_ROOT);
 
%% =========================
%  DT POR MALLA
% ==========================
mesh_list = SETTINGS.MESH_LIST(:).';
ncells_nominal = (mesh_list - 2) / 2;
h_nominal = 1 ./ ncells_nominal;
dt_list = SETTINGS.DT_FINEST * (h_nominal / h_nominal(end));
 
fprintf('=== ESTUDIO EXPLÍCITO ===\n');
for k = 1:numel(mesh_list)
    fprintf('Malla %d -> ncells ~ %d, h ~ %.4e, dt = %.4e\n', ...
        mesh_list(k), ncells_nominal(k), h_nominal(k), dt_list(k));
end
 
%% =========================
%  CASOS EXPLÍCITOS: RUSANOV Y HLLC
% ==========================
rusanov_cases = cell(numel(mesh_list), 1);
hllc_cases = cell(numel(mesh_list), 1);
 
for k = 1:numel(mesh_list)
    fprintf('\n[EXPLÍCITO][RUSANOV] malla %d...\n', mesh_list(k));
    rusanov_cases{k} = simulate_case(mesh_list(k), dt_list(k), 'rusanov', 'explicit', SETTINGS);
 
    fprintf('[EXPLÍCITO][HLLC]    malla %d...\n', mesh_list(k));
    hllc_cases{k} = simulate_case(mesh_list(k), dt_list(k), 'hllc', 'explicit', SETTINGS);
end
 
%% =========================
%  CASOS IMPLÍCITOS OPCIONALES
% ==========================
implicit_cases = [];
if SETTINGS.RUN_IMPLICIT
    fprintf('\n=== ESTUDIO IMPLÍCITO RUSANOV (opcional) ===\n');
    implicit_cases = cell(numel(SETTINGS.IMPLICIT_MESH_LIST), 1);
    for k = 1:numel(SETTINGS.IMPLICIT_MESH_LIST)
        fprintf('[IMPLÍCITO][RUSANOV-FD] malla %d...\n', SETTINGS.IMPLICIT_MESH_LIST(k));
        implicit_cases{k} = simulate_case( ...
            SETTINGS.IMPLICIT_MESH_LIST(k), ...
            SETTINGS.IMPLICIT_DT(k), ...
            'rusanov', ...
            'implicit_fd', ...
            SETTINGS);
    end
end
 
%% =========================
%  SOLUCIÓN EXACTA DE REFERENCIA
% ==========================
verification_case = get_case_by_mesh(hllc_cases, SETTINGS.VERIFICATION_MESH);
x_exact = linspace(0, 1, max(4000, 4 * verification_case.ncells)).';
exact = exact_sod_solution( ...
    x_exact, SETTINGS.T_END, SETTINGS.X0, SETTINGS.LEFT, SETTINGS.RIGHT, SETTINGS.GAMMA);
 
%% =========================
%  FIGURA 1: VERIFICACIÓN HLLC vs EXACTA
% ==========================
fig = figure('Name', 'verification_hllc_exact', 'Color', 'w');
plot_four_variables_against_exact(verification_case, exact, ...
    sprintf('Verificación Sod: HLLC (%d nodos) vs solución exacta', SETTINGS.VERIFICATION_MESH));
save_figure(fig, SETTINGS, 'fig01_verificacion_hllc_vs_exacta');
 
%% =========================
%  FIGURA 2: RUSANOV vs HLLC vs EXACTA
% ==========================
case_rus_cmp = get_case_by_mesh(rusanov_cases, SETTINGS.COMPARISON_MESH);
case_hllc_cmp = get_case_by_mesh(hllc_cases, SETTINGS.COMPARISON_MESH);
fig = figure('Name', 'comparison_rusanov_hllc_exact', 'Color', 'w');
plot_method_comparison(case_rus_cmp, case_hllc_cmp, exact, ...
    sprintf('Comparación de métodos en malla de %d nodos', SETTINGS.COMPARISON_MESH));
save_figure(fig, SETTINGS, 'fig02_rusanov_vs_hllc_vs_exacta');
 
%% =========================
%  FIGURA 3: REFINAMIENTO DE MALLA (HLLC)
% ==========================
fig = figure('Name', 'mesh_refinement_hllc', 'Color', 'w');
plot_mesh_refinement(hllc_cases, mesh_list, SETTINGS.REFINEMENT_PLOT_SKIP, ...
    'Refinamiento de malla con HLLC');
save_figure(fig, SETTINGS, 'fig03_refinamiento_malla_hllc');
 
%% =========================
%  ERRORES L1 FRENTE A LA SOLUCIÓN EXACTA
% ==========================
err_rus = compute_error_struct(rusanov_cases, exact);
err_hllc = compute_error_struct(hllc_cases, exact);
 
% Magnitud suave para GCI/AR por tríos consecutivos, igual que en el script de referencia
Ek_rus = compute_mean_specific_kinetic_energy(rusanov_cases);
Ek_hllc = compute_mean_specific_kinetic_energy(hllc_cases);
 
%% =========================
%  DATOS IMPLÍCITOS PARA FIGURAS 4 Y 5
% ==========================
if SETTINGS.RUN_IMPLICIT && ~isempty(implicit_cases)
    imp_ncells  = cellfun(@(c) c.ncells,  implicit_cases);
    imp_h       = cellfun(@(c) c.h,       implicit_cases);
    imp_runtime = cellfun(@(c) c.runtime, implicit_cases);
    imp_nnodes  = cellfun(@(c) c.nnodes,  implicit_cases);
    imp_nsteps  = cellfun(@(c) c.nsteps_est, implicit_cases);
    imp_dt      = cellfun(@(c) c.dt,      implicit_cases);
 
    % Orden del coste implícito (sobre las mallas disponibles)
    mask_imp = imp_ncells >= SETTINGS.COST_ORDER_MIN_CELLS;
    if sum(mask_imp) >= 2
        cost_ord_imp = compute_cost_order_filtered(imp_ncells(:), imp_h(:), imp_runtime(:), SETTINGS.COST_ORDER_MIN_CELLS);
    else
        cost_ord_imp = [];
    end
    fit_imp = fit_cost_slope(imp_ncells(:), imp_runtime(:), min(imp_ncells));
else
    imp_ncells  = [];
    imp_h       = [];
    imp_runtime = [];
    imp_nnodes  = [];
    imp_nsteps  = [];
    imp_dt      = [];
    cost_ord_imp = [];
    fit_imp = [];
end
 
%% =========================
%  FIGURA 4: COSTE COMPUTACIONAL (EXPLÍCITO + IMPLÍCITO OPCIONAL)
% ==========================
 
% --- Tabla en Command Window ---
fprintf('\n=== TABLA DE COSTES COMPUTACIONALES ===\n');
nMesh_exp = numel(err_rus.ncells);
cost_table_rus = table( ...
    repmat("Rusanov (explícito)", nMesh_exp, 1), ...
    err_rus.ncells(:), err_rus.h(:), ...
    err_rus.dt(:), ...
    err_rus.nsteps(:), err_rus.runtime(:), ...
    'VariableNames', {'Método','Ncells','h','dt_usado','Pasos_est','Runtime_s'});
disp(cost_table_rus);
 
cost_table_hllc = table( ...
    repmat("HLLC (explícito)", nMesh_exp, 1), ...
    err_hllc.ncells(:), err_hllc.h(:), ...
    err_hllc.dt(:), ...
    err_hllc.nsteps(:), err_hllc.runtime(:), ...
    'VariableNames', {'Método','Ncells','h','dt_usado','Pasos_est','Runtime_s'});
disp(cost_table_hllc);
 
if SETTINGS.RUN_IMPLICIT && ~isempty(imp_ncells)
    nMesh_imp = numel(imp_ncells);
    cost_table_imp = table( ...
        repmat("Rusanov (implícito FD)", nMesh_imp, 1), ...
        imp_ncells(:), imp_h(:), imp_dt(:), ...
        imp_nsteps(:), imp_runtime(:), ...
        'VariableNames', {'Método','Ncells','h','dt_usado','Pasos_est','Runtime_s'});
    disp(cost_table_imp);
end
 
% --- Figura ---
fig = figure('Name', 'cost_explicit', 'Color', 'w');
hold on;
 
% Explícitos
plot(err_rus.ncells,  err_rus.runtime,  'o-', 'LineWidth', 2, 'MarkerSize', 8);
plot(err_hllc.ncells, err_hllc.runtime, 's-', 'LineWidth', 2, 'MarkerSize', 8);
 
% Implícito (si existe)
if SETTINGS.RUN_IMPLICIT && ~isempty(imp_ncells)
    plot(imp_ncells, imp_runtime, 'd-', 'LineWidth', 2, 'MarkerSize', 8);
end
 
% Referencias cuadráticas ancladas a la primera malla >= COST_ORDER_MIN_CELLS
idx_anchor_rus  = find(err_rus.ncells  >= SETTINGS.COST_ORDER_MIN_CELLS, 1, 'first');
idx_anchor_hllc = find(err_hllc.ncells >= SETTINGS.COST_ORDER_MIN_CELLS, 1, 'first');
nref = linspace(min(err_rus.ncells), max(err_rus.ncells), 300).';
ref_rus  = err_rus.runtime(idx_anchor_rus)  * (nref / err_rus.ncells(idx_anchor_rus)).^2;
ref_hllc = err_hllc.runtime(idx_anchor_hllc) * (nref / err_hllc.ncells(idx_anchor_hllc)).^2;
plot(nref, ref_rus,  '--', 'LineWidth', 1.2, 'Color', [0.3 0.3 0.9]);
plot(nref, ref_hllc, '--', 'LineWidth', 1.2, 'Color', [0.9 0.3 0.3]);
 
grid on;
xlabel('Número de celdas [-]');
ylabel('Tiempo de cálculo [s]');
title('Coste computacional');
 
if SETTINGS.RUN_IMPLICIT && ~isempty(imp_ncells)
    legend('Rusanov (explícito)', 'HLLC (explícito)', 'Rusanov (implícito FD)', ...
        'Ref. cuadrática (Rusanov exp.)', 'Ref. cuadrática (HLLC exp.)', ...
        'Location', 'northwest');
else
    legend('Rusanov', 'HLLC', ...
        'Ref. cuadrática (Rusanov)', 'Ref. cuadrática (HLLC)', ...
        'Location', 'northwest');
end
save_figure(fig, SETTINGS, 'fig04_coste_computacional');
 
%% =========================
%  FIGURA 5: ORDEN DEL COSTE COMPUTACIONAL (EXPLÍCITO + IMPLÍCITO OPCIONAL)
% ==========================
cost_ord_rus  = compute_cost_order_filtered(err_rus.ncells,  err_rus.h,  err_rus.runtime,  SETTINGS.COST_ORDER_MIN_CELLS);
cost_ord_hllc = compute_cost_order_filtered(err_hllc.ncells, err_hllc.h, err_hllc.runtime, SETTINGS.COST_ORDER_MIN_CELLS);
fit_rus  = fit_cost_slope(err_rus.ncells,  err_rus.runtime,  SETTINGS.COST_ORDER_MIN_CELLS);
fit_hllc = fit_cost_slope(err_hllc.ncells, err_hllc.runtime, SETTINGS.COST_ORDER_MIN_CELLS);
 
% --- Tabla en Command Window ---
fprintf('\n=== TABLA DE ORDEN DEL COSTE COMPUTACIONAL ===\n');
order_table_rus = table( ...
    cost_ord_rus.h_mid(:), cost_ord_rus.p(:), ...
    'VariableNames', {'h_mid_Rusanov', 'Orden_Rusanov'});
order_table_hllc = table( ...
    cost_ord_hllc.h_mid(:), cost_ord_hllc.p_plot(:), ...
    'VariableNames', {'h_mid_HLLC', 'Orden_HLLC'});
fprintf('  Rusanov (ajuste global p = %.4f):\n', fit_rus.p);
disp(order_table_rus);
fprintf('  HLLC    (ajuste global p = %.4f):\n', fit_hllc.p);
disp(order_table_hllc);
 
if SETTINGS.RUN_IMPLICIT && ~isempty(cost_ord_imp)
    fprintf('  Rusanov implícito (ajuste global p = %.4f):\n', fit_imp.p);
    order_table_imp = table( ...
        cost_ord_imp.h_mid(:), cost_ord_imp.p_plot(:), ...
        'VariableNames', {'h_mid_Implícito', 'Orden_Implícito'});
    disp(order_table_imp);
end
 
% --- Figura ---
if SETTINGS.RUN_IMPLICIT && ~isempty(cost_ord_imp)
    title_str = sprintf('Orden del coste (ajuste: Rusanov %.2f, HLLC %.2f, Impl. %.2f)', ...
        fit_rus.p, fit_hllc.p, fit_imp.p);
else
    title_str = sprintf('Orden del coste (ajuste global: Rusanov %.2f, HLLC %.2f)', ...
        fit_rus.p, fit_hllc.p);
end
 
fig = figure('Name', 'cost_order', 'Color', 'w');
hold on;
semilogx(cost_ord_rus.h_plot,  cost_ord_rus.p_plot,  'o-', 'LineWidth', 2, 'MarkerSize', 8);
semilogx(cost_ord_hllc.h_plot, cost_ord_hllc.p_plot, 's-', 'LineWidth', 2, 'MarkerSize', 8);
if SETTINGS.RUN_IMPLICIT && ~isempty(cost_ord_imp)
    semilogx(cost_ord_imp.h_plot, cost_ord_imp.p_plot, 'd-', 'LineWidth', 2, 'MarkerSize', 8);
end
yline(2, '--k', 'LineWidth', 1.2);
grid on;
xlabel('h [-]');
ylabel('Orden del coste [-]');
title(title_str);
if SETTINGS.RUN_IMPLICIT && ~isempty(cost_ord_imp)
    legend('Rusanov (explícito)', 'HLLC (explícito)', 'Rusanov (implícito FD)', ...
        'Orden 2 (referencia)', 'Location', 'best');
else
    legend('Rusanov', 'HLLC', 'Orden 2 (referencia)', 'Location', 'best');
end
save_figure(fig, SETTINGS, 'fig05_orden_coste_computacional');
 
%% =========================
%  FIGURA 6: ERROR L1 vs h
% ==========================
fig = figure('Name', 'error_vs_h', 'Color', 'w');
tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact');
vars = {'rho','u','p'};
ylabels = {'Error L_1(\rho)', 'Error L_1(u)', 'Error L_1(p)'};
for i = 1:3
    nexttile;
    loglog(err_rus.h, err_rus.(vars{i}), 'o-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
    loglog(err_hllc.h, err_hllc.(vars{i}), 's-', 'LineWidth', 1.8, 'MarkerSize', 7);
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('h [-]');
    ylabel(ylabels{i});
    title(['Convergencia espacial de ', vars{i}]);
    legend('Rusanov', 'HLLC', 'Location', 'best');
end
save_figure(fig, SETTINGS, 'fig06_error_L1_vs_h');
 
%% =========================
%  FIGURA 7: ORDEN OBSERVADO
% ==========================
ord_rus = compute_observed_order(err_rus.h, err_rus);
ord_hllc = compute_observed_order(err_hllc.h, err_hllc);
 
fig = figure('Name', 'observed_order', 'Color', 'w');
tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:3
    nexttile;
    semilogx(ord_rus.h_mid, ord_rus.(vars{i}), 'o-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
    semilogx(ord_hllc.h_mid, ord_hllc.(vars{i}), 's-', 'LineWidth', 1.8, 'MarkerSize', 7);
    yline(SETTINGS.ORDER_REFERENCE, '--', 'LineWidth', 1.1);
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabel('h medio [-]');
    ylabel('p observado [-]');
    title(['Orden observado en ', vars{i}]);
    legend('Rusanov', 'HLLC', sprintf('Referencia p=%.0f', SETTINGS.ORDER_REFERENCE), 'Location', 'best');
end
save_figure(fig, SETTINGS, 'fig07_orden_observado');
 
%% =========================
%  FIGURAS 8 y 9: RICHARDSON + EXTRAPOLACIÓN EN DENSIDAD
% ==========================
rich_rus = compute_richardson_on_three_finest(rusanov_cases);
rich_hllc = compute_richardson_on_three_finest(hllc_cases);
 
fig = figure('Name', 'richardson_rusanov_density', 'Color', 'w');
plot_richardson_density(rich_rus, 'Richardson en densidad - Rusanov');
save_figure(fig, SETTINGS, 'fig08_richardson_densidad_rusanov');
 
fig = figure('Name', 'richardson_hllc_density', 'Color', 'w');
plot_richardson_density(rich_hllc, 'Richardson en densidad - HLLC');
save_figure(fig, SETTINGS, 'fig09_richardson_densidad_hllc');
 
%% =========================
%  FIGURAS 10-13: GCI / AR POR TRÍOS CONSECUTIVOS SOBRE ENERGÍA CINÉTICA
%  ESPECÍFICA MEDIA, IGUAL QUE EN EL SCRIPT DE REFERENCIA
% ==========================
gci_series_rus = compute_gci_series_from_scalar(err_rus.h, Ek_rus);
gci_series_hllc = compute_gci_series_from_scalar(err_hllc.h, Ek_hllc);
 
fig = figure('Name', 'ar_rusanov', 'Color', 'w');
plot_ar_series_single(gci_series_rus, 'Rango asintótico (AR) por trío de mallas - Rusanov');
save_figure(fig, SETTINGS, 'fig10_ar_rusanov');
 
fig = figure('Name', 'ar_hllc', 'Color', 'w');
plot_ar_series_single(gci_series_hllc, 'Rango asintótico (AR) por trío de mallas - HLLC');
save_figure(fig, SETTINGS, 'fig11_ar_hllc');
 
fig = figure('Name', 'gci_rusanov', 'Color', 'w');
plot_gci_series_single(gci_series_rus, 'Índice de convergencia de malla - Rusanov');
save_figure(fig, SETTINGS, 'fig12_gci_rusanov');
 
fig = figure('Name', 'gci_hllc', 'Color', 'w');
plot_gci_series_single(gci_series_hllc, 'Índice de convergencia de malla - HLLC');
save_figure(fig, SETTINGS, 'fig13_gci_hllc');
 
 
%% =========================
%  TABLAS EN .MAT
% ==========================
if SETTINGS.SAVE_MAT
    save(fullfile(SETTINGS.SAVE_DIR, 'report_figure_data_v2.mat'), ...
        'SETTINGS', 'rusanov_cases', 'hllc_cases', 'implicit_cases', ...
        'exact', 'err_rus', 'err_hllc', 'ord_rus', 'ord_hllc', ...
        'rich_rus', 'rich_hllc', 'gci_series_rus', 'gci_series_hllc', 'cost_ord_rus', 'cost_ord_hllc', 'fit_rus', 'fit_hllc');
end
 
fprintf('\nFiguras guardadas en:\n%s\n', SETTINGS.SAVE_DIR);
 
%% ========================================================================
%  FUNCIONES LOCALES
% ========================================================================
function setup_project_paths(project_root)
    addpath(fullfile(project_root, 'mesh_processing'));
    addpath(genpath(fullfile(project_root, 'utils')));
    addpath(fullfile(project_root, 'constants'));
    addpath(fullfile(project_root, 'initial_conditions'));
    addpath(fullfile(project_root, 'problems', 'fvm_1D_euler'));
    addpath(genpath(fullfile(project_root, 'convective_flux')));
    addpath(fullfile(project_root, 'stopping_criteria'));
    addpath(fullfile(project_root, 'timestep_control'));
    addpath(fullfile(project_root, 'propagators'));
    addpath(fullfile(project_root, 'results_manager'));
    addpath(fullfile(project_root, 'sources'));
    addpath(fullfile(project_root, 'boundary_condition'));
    addpath(fullfile(project_root, 'convergency_study'));
end
 
function out = simulate_case(nnodes, dt_value, flux_name, time_mode, SETTINGS)
    [nodes_file, cells_file, bc_files] = mesh_case_paths(nnodes);
    cells = mesh_processor(nodes_file, cells_file, bc_files);

    centroids_x = get_cell_centroid_x(cells);
    centroids_x = centroids_x(:);
    num_cells = numel(centroids_x);

    w0 = Config.INITIAL_CONDITIONS(centroids_x);

    switch lower(time_mode)
        case 'explicit'
            propagator = @fw_euler_sparse;
            problem = @(state, time) local_fvm_problem(state, cells, time, flux_name, SETTINGS.GAMMA);
        case 'implicit_fd'
            propagator = @bw_euler_sparse;
            problem = @(state, time) fvm_1D_euler_implicit(state, cells);
        otherwise
            error('Modo temporal no reconocido: %s', time_mode);
    end

    runtime_samples = zeros(SETTINGS.RUNTIME_REPEATS, 1);
    nsteps_samples  = zeros(SETTINGS.RUNTIME_REPEATS, 1);

    for irun = 1:SETTINGS.RUNTIME_REPEATS
        w_run = w0;
        t_run = SETTINGS.T0;
        nsteps_run = 0;

        tic;
        while t_run < SETTINGS.T_END - 1e-14
            dt = min(dt_value, SETTINGS.T_END - t_run);
            w_run = propagator(w_run, t_run, dt, problem);
            t_run = t_run + dt;
            nsteps_run = nsteps_run + 1;
        end
        runtime_samples(irun) = toc;
        nsteps_samples(irun) = nsteps_run;
    end

    if SETTINGS.RUNTIME_DROP_FIRST && numel(runtime_samples) >= 3
        runtime = median(runtime_samples(2:end));
    else
        runtime = median(runtime_samples);
    end

    nsteps_est = round(median(nsteps_samples));

    w = w_run;
    t = t_run;

    rho  = w(1:num_cells);
    rhou = w(num_cells+1:2*num_cells);
    E    = w(2*num_cells+1:3*num_cells);

    u = rhou ./ rho;
    p = (SETTINGS.GAMMA - 1) .* (E - 0.5 .* (rhou.^2) ./ rho);

    [x_sorted, idx] = sort(centroids_x);
    rho = rho(idx);
    u   = u(idx);
    p   = p(idx);
    E   = E(idx);
    ei  = p ./ ((SETTINGS.GAMMA - 1) .* rho);

    out = struct();
    out.nnodes = nnodes;
    out.ncells = num_cells;
    out.h = mean(diff(x_sorted));
    out.dt = dt_value;
    out.t = t;
    out.runtime = runtime;
    out.runtime_samples = runtime_samples;
    out.nsteps_est = nsteps_est;
    out.nsteps_samples = nsteps_samples;
    out.x = x_sorted;
    out.rho = rho;
    out.u = u;
    out.p = p;
    out.E = E;
    out.ei = ei;
    out.flux = flux_name;
    out.time_mode = time_mode;
end
 
function [A, b] = local_fvm_problem(state, cells, ~, flux_name, gamma)
    num_cells = numel(cells);
    A = sparse(3*num_cells, 3*num_cells);
    b = local_residual_1d(state, cells, flux_name, gamma);
end
 
function b = local_residual_1d(state, cells, flux_name, gamma)
    num_cells = numel(cells);
    [x_sorted, perm, rho, rhou, E] = sort_state_by_cell_centroid_x(state, cells);
 
    rhs_rho = zeros(num_cells, 1);
    rhs_rhou = zeros(num_cells, 1);
    rhs_E = zeros(num_cells, 1);
 
    for i = 1:num_cells
        Ui = [rho(i); rhou(i); E(i)];
 
        if i == 1
            UL = closed_bc(Ui);
            xL = x_sorted(i) - 0.5 * (x_sorted(i+1) - x_sorted(i));
        else
            UL = [rho(i-1); rhou(i-1); E(i-1)];
            xL = 0.5 * (x_sorted(i-1) + x_sorted(i));
        end
 
        if i == num_cells
            UR = closed_bc(Ui);
            xR = x_sorted(i) + 0.5 * (x_sorted(i) - x_sorted(i-1));
        else
            UR = [rho(i+1); rhou(i+1); E(i+1)];
            xR = 0.5 * (x_sorted(i) + x_sorted(i+1));
        end
 
        switch lower(flux_name)
            case 'rusanov'
                F_left = rusanov_flux_local(UL, Ui, gamma);
                F_right = rusanov_flux_local(Ui, UR, gamma);
            case 'hllc'
                F_left = hllc_flux(UL, Ui, gamma);
                F_right = hllc_flux(Ui, UR, gamma);
            otherwise
                error('Flujo no reconocido: %s', flux_name);
        end
 
        dx_i = xR - xL;
        rhs_i = -(F_right - F_left) / dx_i;
 
        rhs_rho(i) = rhs_i(1);
        rhs_rhou(i) = rhs_i(2);
        rhs_E(i) = rhs_i(3);
    end
 
    rhs_rho_orig = zeros(num_cells, 1);
    rhs_rhou_orig = zeros(num_cells, 1);
    rhs_E_orig = zeros(num_cells, 1);
 
    rhs_rho_orig(perm) = rhs_rho;
    rhs_rhou_orig(perm) = rhs_rhou;
    rhs_E_orig(perm) = rhs_E;
 
    b = [rhs_rho_orig; rhs_rhou_orig; rhs_E_orig];
end
 
function F = rusanov_flux_local(UL, UR, gamma)
    FL = euler_flux_1d(UL, gamma);
    FR = euler_flux_1d(UR, gamma);
 
    aL = max_wave_speed_local(UL, gamma);
    aR = max_wave_speed_local(UR, gamma);
    alpha = max(aL, aR);
 
    F = 0.5 * (FL + FR) - 0.5 * alpha * (UR - UL);
end
 
function a = max_wave_speed_local(U, gamma)
    [rho, u, p, ~] = cons2prim(U, gamma);
    c = sqrt(gamma * p / rho);
    a = abs(u) + c;
end
 
function case_out = get_case_by_mesh(cases_in, mesh_tag)
    idx = find(cellfun(@(c) c.nnodes == mesh_tag, cases_in), 1, 'first');
    if isempty(idx)
        error('No se encontró la malla %d.', mesh_tag);
    end
    case_out = cases_in{idx};
end
 
function plot_four_variables_against_exact(case_num, exact, fig_title)
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
 
    nexttile;
    plot(exact.x, exact.rho, 'k-', 'LineWidth', 1.4); hold on;
    plot(case_num.x, case_num.rho, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.8);
    grid on; xlabel('x [-]'); ylabel('\rho [-]'); title('Densidad');
    legend('Exacta', 'Numérica', 'Location', 'best');
 
    nexttile;
    plot(exact.x, exact.u, 'k-', 'LineWidth', 1.4); hold on;
    plot(case_num.x, case_num.u, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.8);
    grid on; xlabel('x [-]'); ylabel('u [-]'); title('Velocidad');
    legend('Exacta', 'Numérica', 'Location', 'best');
 
    nexttile;
    plot(exact.x, exact.p, 'k-', 'LineWidth', 1.4); hold on;
    plot(case_num.x, case_num.p, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.8);
    grid on; xlabel('x [-]'); ylabel('p [-]'); title('Presión');
    legend('Exacta', 'Numérica', 'Location', 'best');
 
    nexttile;
    plot(exact.x, exact.ei, 'k-', 'LineWidth', 1.4); hold on;
    plot(case_num.x, case_num.ei, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.8);
    grid on; xlabel('x [-]'); ylabel('e_i [-]'); title('Energía interna específica');
    legend('Exacta', 'Numérica', 'Location', 'best');
 
    sgtitle(fig_title, 'FontWeight', 'bold');
end
 
function plot_method_comparison(case_rus, case_hllc, exact, fig_title)
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');
 
    nexttile;
    plot(exact.x, exact.rho, 'k-', 'LineWidth', 1.4); hold on;
    stairs(case_rus.x, case_rus.rho, '-', 'LineWidth', 1.2);
    stairs(case_hllc.x, case_hllc.rho, '-', 'LineWidth', 1.2);
    grid on; xlabel('x [-]'); ylabel('\rho [-]'); title('Densidad');
    legend('Exacta', 'Rusanov', 'HLLC', 'Location', 'best');
 
    nexttile;
    plot(exact.x, exact.u, 'k-', 'LineWidth', 1.4); hold on;
    stairs(case_rus.x, case_rus.u, '-', 'LineWidth', 1.2);
    stairs(case_hllc.x, case_hllc.u, '-', 'LineWidth', 1.2);
    grid on; xlabel('x [-]'); ylabel('u [-]'); title('Velocidad');
    legend('Exacta', 'Rusanov', 'HLLC', 'Location', 'best');
 
    nexttile;
    plot(exact.x, exact.p, 'k-', 'LineWidth', 1.4); hold on;
    stairs(case_rus.x, case_rus.p, '-', 'LineWidth', 1.2);
    stairs(case_hllc.x, case_hllc.p, '-', 'LineWidth', 1.2);
    grid on; xlabel('x [-]'); ylabel('p [-]'); title('Presión');
    legend('Exacta', 'Rusanov', 'HLLC', 'Location', 'best');
 
    nexttile;
    plot(exact.x, exact.ei, 'k-', 'LineWidth', 1.4); hold on;
    stairs(case_rus.x, case_rus.ei, '-', 'LineWidth', 1.2);
    stairs(case_hllc.x, case_hllc.ei, '-', 'LineWidth', 1.2);
    grid on; xlabel('x [-]'); ylabel('e_i [-]'); title('Energía interna específica');
    legend('Exacta', 'Rusanov', 'HLLC', 'Location', 'best');
 
    sgtitle(fig_title, 'FontWeight', 'bold');
end
 
function plot_mesh_refinement(cases_in, mesh_list, skip_list, fig_title)
    idx_plot = ~ismember(mesh_list, skip_list);
 
    tiledlayout(3,1, 'Padding', 'compact', 'TileSpacing', 'compact');
 
    nexttile; hold on;
    for k = find(idx_plot)
        stairs(cases_in{k}.x, cases_in{k}.rho, 'LineWidth', 1.4);
    end
    grid on; xlabel('x [m]'); ylabel('\rho [kg/m^3]'); title('Densidad');
    legend(compose('%d nodos', mesh_list(idx_plot)), 'Location', 'best');
 
    nexttile; hold on;
    for k = find(idx_plot)
        stairs(cases_in{k}.x, cases_in{k}.u, 'LineWidth', 1.4);
    end
    grid on; xlabel('x [m]'); ylabel('u [m/s]'); title('Velocidad');
    legend(compose('%d nodos', mesh_list(idx_plot)), 'Location', 'best');
 
    nexttile; hold on;
    for k = find(idx_plot)
        stairs(cases_in{k}.x, cases_in{k}.p, 'LineWidth', 1.4);
    end
    grid on; xlabel('x [m]'); ylabel('p [Pa]'); title('Presión');
    legend(compose('%d nodos', mesh_list(idx_plot)), 'Location', 'best');
 
    sgtitle(fig_title, 'FontWeight', 'bold');
end
 
function err = compute_error_struct(cases_in, exact)
    n = numel(cases_in);

    err.h = zeros(n,1);
    err.ncells = zeros(n,1);
    err.nnodes = zeros(n,1);
    err.dt = zeros(n,1);
    err.nsteps = zeros(n,1);
    err.runtime = zeros(n,1);
    err.rho = zeros(n,1);
    err.u = zeros(n,1);
    err.p = zeros(n,1);

    for k = 1:n
        x = cases_in{k}.x;
        rho_ex = interp1(exact.x, exact.rho, x, 'linear');
        u_ex   = interp1(exact.x, exact.u, x, 'linear');
        p_ex   = interp1(exact.x, exact.p, x, 'linear');

        err.h(k) = cases_in{k}.h;
        err.ncells(k) = cases_in{k}.ncells;
        err.nnodes(k) = cases_in{k}.nnodes;
        err.dt(k) = cases_in{k}.dt;
        err.nsteps(k) = cases_in{k}.nsteps_est;
        err.runtime(k) = cases_in{k}.runtime;

        err.rho(k) = mean(abs(cases_in{k}.rho - rho_ex));
        err.u(k)   = mean(abs(cases_in{k}.u   - u_ex));
        err.p(k)   = mean(abs(cases_in{k}.p   - p_ex));
    end
end
 
function ord = compute_observed_order(h, err)
    ord.h_mid = sqrt(h(1:end-1) .* h(2:end));
    ord.rho = abs(log(err.rho(2:end) ./ err.rho(1:end-1)) ./ log(h(2:end) ./ h(1:end-1)));
    ord.u   = abs(log(err.u(2:end) ./ err.u(1:end-1)) ./ log(h(2:end) ./ h(1:end-1)));
    ord.p   = abs(log(err.p(2:end) ./ err.p(1:end-1)) ./ log(h(2:end) ./ h(1:end-1)));
end
 
function cost_ord = compute_cost_order_filtered(ncells, h, runtime, min_cells)
    mask = ncells >= min_cells;
    N = ncells(mask);
    h_use = h(mask);
    T = runtime(mask);
 
    p = abs(log(T(2:end) ./ T(1:end-1)) ./ log(h_use(2:end) ./ h_use(1:end-1)));
    h_mid = sqrt(h_use(1:end-1) .* h_use(2:end));
 
    cost_ord.N = N(:);
    cost_ord.h = h_use(:);
    cost_ord.h_mid = h_mid(:);
    cost_ord.p = p(:);
    cost_ord.h_plot = flip(h_mid(:));
    cost_ord.p_plot = flip(p(:));
end
 
function fit_out = fit_cost_slope(ncells, runtime, min_cells)
    mask = ncells >= min_cells;
    N = ncells(mask);
    T = runtime(mask);
 
    coeffs = polyfit(log(N), log(T), 1);
    fit_out.p = coeffs(1);
    fit_out.a = coeffs(2);
    fit_out.N = N(:);
    fit_out.Tfit = exp(polyval(coeffs, log(N(:))));
end
 
function rich = compute_richardson_on_three_finest(cases_in)
    coarse = cases_in{end-2};
    medium = cases_in{end-1};
    fine   = cases_in{end};
 
    rho_metrics = richardson_gci_three_grids( ...
        coarse.x, coarse.rho, medium.x, medium.rho, fine.x, fine.rho);
 
    rich = struct();
    rich.coarse = coarse;
    rich.medium = medium;
    rich.fine = fine;
    rich.rho = rho_metrics;
end
 
function Ek_mean = compute_mean_specific_kinetic_energy(cases_in)
    n = numel(cases_in);
    Ek_mean = zeros(n,1);
    for k = 1:n
        u_k = cases_in{k}.u(:);
        Ek_mean(k) = mean(0.5 * u_k.^2);
    end
end
 
function series = compute_gci_series_from_scalar(h, phi)
    Fs = 1.25;
    h = h(:);
    phi = phi(:);
 
    if numel(h) ~= numel(phi)
        error('h y phi deben tener la misma longitud.');
    end
    if numel(h) < 3
        error('Se necesitan al menos tres mallas para calcular GCI/AR.');
    end
    if any(diff(h) >= 0)
        error('Los h deben ir en orden estrictamente decreciente.');
    end
 
    nTrios = numel(h) - 2;
    series.triplet_index = (1:nTrios).';
    series.label_h = strings(nTrios,1);
    series.GCI32 = nan(nTrios,1);
    series.GCI21 = nan(nTrios,1);
    series.AR = nan(nTrios,1);
 
    for i = 1:nTrios
        h3 = h(i);   h2 = h(i+1); h1 = h(i+2);
        phi3 = phi(i); phi2 = phi(i+1); phi1 = phi(i+2);
 
        r32 = h3 / h2;
        r21 = h2 / h1;
        e32 = phi3 - phi2;
        e21 = phi2 - phi1;
        ratio = e32 / e21;
 
        if ratio > 0
            s = 1;
        elseif ratio < 0
            s = -1;
        else
            s = 0;
        end
 
        if abs(e21) < eps || abs(e32) < eps
            p = NaN;
        else
            p0 = abs(log(abs(e32/e21)) / log(r21));
            if ~isfinite(p0) || p0 == 0
                p0 = 1;
            end
            p = p0;
            for kiter = 1:100
                num = r21^p - s;
                den = r32^p - s;
                if abs(den) < eps || abs(num) < eps
                    p = NaN;
                    break;
                end
                p_new = abs((log(abs(e32/e21)) + log(abs(num/den))) / log(r21));
                if ~isfinite(p_new)
                    p = NaN;
                    break;
                end
                if abs(p_new - p) < 1e-10
                    p = p_new;
                    break;
                end
                p = p_new;
            end
        end
 
        ea32 = abs((phi2 - phi3) / phi2);
        ea21 = abs((phi1 - phi2) / phi1);
 
        if isnan(p) || abs(r32^p - 1) < eps || abs(r21^p - 1) < eps
            GCI32 = NaN;
            GCI21 = NaN;
            AR = NaN;
        else
            GCI32 = Fs * ea32 / (r32^p - 1);
            GCI21 = Fs * ea21 / (r21^p - 1);
            if abs(GCI32) < eps
                AR = NaN;
            else
                AR = (GCI21 * r21^p) / GCI32;
            end
        end
 
        series.label_h(i) = sprintf('[%.4g, %.4g, %.4g]', h3, h2, h1);
        series.GCI32(i) = GCI32;
        series.GCI21(i) = GCI21;
        series.AR(i) = AR;
    end
end
 
function plot_richardson_density(rich, fig_title)
    stairs(rich.coarse.x, rich.rho.phi3, 'LineWidth', 1.4); hold on;
    stairs(rich.coarse.x, rich.rho.phi2_on_3, 'LineWidth', 1.4);
    stairs(rich.coarse.x, rich.rho.phi1_on_3, 'LineWidth', 1.4);
    plot(rich.coarse.x, rich.rho.phi_ext, '--', 'LineWidth', 1.8);
    grid on;
    xlabel('x [m]'); ylabel('\rho [kg/m^3]');
    title(fig_title);
    legend( ...
        sprintf('%d nodos', rich.coarse.nnodes), ...
        sprintf('%d nodos \rightarrow malla gruesa', rich.medium.nnodes), ...
        sprintf('%d nodos \rightarrow malla gruesa', rich.fine.nnodes), ...
        'Extrapolación de Richardson', ...
        'Location', 'best');
end
 
function plot_ar_series_single(series, fig_title)
    x = 1:numel(series.AR);
    plot(x, series.AR, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
    yline(1, '--', 'LineWidth', 1.2);
    grid on;
    xlabel('Índice del trío de mallas');
    ylabel('AR');
    title(fig_title);
    xticks(x);
    xticklabels(cellstr(series.label_h));
    xtickangle(20);
end
 
function plot_gci_series_single(series, fig_title)
    x = 1:numel(series.GCI32);
    plot(x, 100*series.GCI32, 's-', 'LineWidth', 1.8, 'MarkerSize', 7); hold on;
    plot(x, 100*series.GCI21, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
    grid on;
    xlabel('Índice del trío de mallas');
    ylabel('GCI [%]');
    title(fig_title);
    legend('GCI_{32}', 'GCI_{21}', 'Location', 'best');
    xticks(x);
    xticklabels(cellstr(series.label_h));
    xtickangle(20);
end
 
function save_figure(fig, SETTINGS, base_name)
    if ~SETTINGS.SAVE_FIGURES
        return;
    end
    png_file = fullfile(SETTINGS.SAVE_DIR, [base_name, '.png']);
    fig_file = fullfile(SETTINGS.SAVE_DIR, [base_name, '.fig']);
    print(fig, png_file, '-dpng', '-r300');
    savefig(fig, fig_file);
end
 
function sol = exact_sod_solution(x, t, x0, left, right, gamma)
    if t <= 0
        sol.x = x(:);
        sol.rho = left.rho * ones(size(x(:)));
        sol.u = left.u * ones(size(x(:)));
        sol.p = left.p * ones(size(x(:)));
        right_mask = x(:) > x0;
        sol.rho(right_mask) = right.rho;
        sol.u(right_mask) = right.u;
        sol.p(right_mask) = right.p;
        sol.E = sol.p/(gamma-1) + 0.5*sol.rho.*sol.u.^2;
        sol.ei = sol.p ./ ((gamma - 1) .* sol.rho);
        return;
    end
 
    [p_star, u_star] = star_region_pressure_velocity(left, right, gamma);
    xi = (x(:) - x0) / t;
 
    rho = zeros(size(xi));
    u = zeros(size(xi));
    p = zeros(size(xi));
 
    aL = sqrt(gamma * left.p / left.rho);
    aR = sqrt(gamma * right.p / right.rho);
 
    if p_star > left.p
        rho_star_L = left.rho * ((p_star/left.p + (gamma-1)/(gamma+1)) / ...
            ((gamma-1)/(gamma+1) * p_star/left.p + 1));
    else
        rho_star_L = left.rho * (p_star/left.p)^(1/gamma);
    end
 
    if p_star > right.p
        rho_star_R = right.rho * ((p_star/right.p + (gamma-1)/(gamma+1)) / ...
            ((gamma-1)/(gamma+1) * p_star/right.p + 1));
    else
        rho_star_R = right.rho * (p_star/right.p)^(1/gamma);
    end
 
    for i = 1:numel(xi)
        s = xi(i);
 
        if s <= u_star
            if p_star > left.p
                SL = left.u - aL * sqrt((gamma+1)/(2*gamma) * (p_star/left.p) + (gamma-1)/(2*gamma));
                if s <= SL
                    rho(i) = left.rho; u(i) = left.u; p(i) = left.p;
                else
                    rho(i) = rho_star_L; u(i) = u_star; p(i) = p_star;
                end
            else
                SHL = left.u - aL;
                a_star_L = aL * (p_star/left.p)^((gamma-1)/(2*gamma));
                STL = u_star - a_star_L;
                if s <= SHL
                    rho(i) = left.rho; u(i) = left.u; p(i) = left.p;
                elseif s >= STL
                    rho(i) = rho_star_L; u(i) = u_star; p(i) = p_star;
                else
                    ufan = 2/(gamma+1) * (aL + 0.5*(gamma-1)*left.u + s);
                    afan = 2/(gamma+1) * (aL + 0.5*(gamma-1)*(left.u - s));
                    rho(i) = left.rho * (afan/aL)^(2/(gamma-1));
                    u(i) = ufan;
                    p(i) = left.p * (afan/aL)^(2*gamma/(gamma-1));
                end
            end
        else
            if p_star > right.p
                SR = right.u + aR * sqrt((gamma+1)/(2*gamma) * (p_star/right.p) + (gamma-1)/(2*gamma));
                if s >= SR
                    rho(i) = right.rho; u(i) = right.u; p(i) = right.p;
                else
                    rho(i) = rho_star_R; u(i) = u_star; p(i) = p_star;
                end
            else
                SHR = right.u + aR;
                a_star_R = aR * (p_star/right.p)^((gamma-1)/(2*gamma));
                STR = u_star + a_star_R;
                if s >= SHR
                    rho(i) = right.rho; u(i) = right.u; p(i) = right.p;
                elseif s <= STR
                    rho(i) = rho_star_R; u(i) = u_star; p(i) = p_star;
                else
                    ufan = 2/(gamma+1) * (-aR + 0.5*(gamma-1)*right.u + s);
                    afan = 2/(gamma+1) * (aR - 0.5*(gamma-1)*(right.u - s));
                    rho(i) = right.rho * (afan/aR)^(2/(gamma-1));
                    u(i) = ufan;
                    p(i) = right.p * (afan/aR)^(2*gamma/(gamma-1));
                end
            end
        end
    end
 
    E = p/(gamma-1) + 0.5 * rho .* u.^2;
 
    sol = struct();
    sol.x = x(:);
    sol.rho = rho;
    sol.u = u;
    sol.p = p;
    sol.E = E;
    sol.ei = p ./ ((gamma - 1) .* rho);
end
 
function [p_star, u_star] = star_region_pressure_velocity(left, right, gamma)
    pL = left.p; pR = right.p;
    uL = left.u; uR = right.u;
    rhoL = left.rho; rhoR = right.rho;
 
    aL = sqrt(gamma * pL / rhoL);
    aR = sqrt(gamma * pR / rhoR);
 
    p_guess = max(1e-8, 0.5*(pL + pR) - 0.125*(uR - uL)*(rhoL + rhoR)*(aL + aR));
    p_old = p_guess;
 
    for iter = 1:100
        [fL, dfL] = pressure_function(p_old, left, gamma);
        [fR, dfR] = pressure_function(p_old, right, gamma);
 
        p_new = p_old - (fL + fR + uR - uL) / (dfL + dfR);
        p_new = max(p_new, 1e-10);
 
        if abs(p_new - p_old) / (0.5 * (p_new + p_old)) < 1e-8
            p_old = p_new;
            break;
        end
        p_old = p_new;
    end
 
    p_star = p_old;
    [fL, ~] = pressure_function(p_star, left, gamma);
    [fR, ~] = pressure_function(p_star, right, gamma);
    u_star = 0.5 * (uL + uR + fR - fL);
end
 
function [f, df] = pressure_function(p, state, gamma)
    pk = state.p;
    rhok = state.rho;
    ak = sqrt(gamma * pk / rhok);
 
    if p > pk
        A = 2 / ((gamma + 1) * rhok);
        B = (gamma - 1) / (gamma + 1) * pk;
        sqrt_term = sqrt(A / (p + B));
        f = (p - pk) * sqrt_term;
        df = sqrt_term * (1 - 0.5 * (p - pk) / (p + B));
    else
        expo = (gamma - 1) / (2 * gamma);
        f = 2 * ak / (gamma - 1) * ((p / pk)^expo - 1);
        df = (1 / (rhok * ak)) * (p / pk)^(-(gamma + 1) / (2 * gamma));
    end
end