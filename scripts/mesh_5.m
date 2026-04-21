%% CREATE_MISSING_REPORT_FIGURES_V7
% Genera las figuras principales para el informe del proyecto de Euler 1D.
% Versión revisada:
% - Coste computacional frente al número de celdas
% - Referencia de orden p=1 para el caso de Sod
% - GCI y rango asintótico calculados sobre TODOS los tripletes consecutivos
% - Gráficas GCI/AR en formato línea, no barras
% - Figura adicional del orden del coste computacional
% - En verificación y comparación se usa energía interna específica
% - Coste computacional y orden del coste calculados al estilo benchmark:
%   * t vs Ncells en log-log
%   * orden local usando mallas consecutivas
%   * ajuste global t ~ C*N^p

clc; clear; close all;

%% =========================
%  SETTINGS BÁSICOS
% ==========================
SETTINGS.PROJECT_ROOT = fileparts(mfilename('fullpath'));

% Paths del proyecto (inline para evitar problemas con setup_project_paths)
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'mesh_processing'));
addpath(genpath(fullfile(SETTINGS.PROJECT_ROOT, 'utils')));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'constants'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'initial_conditions'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'problems', 'fvm_1D_euler'));
addpath(genpath(fullfile(SETTINGS.PROJECT_ROOT, 'convective_flux')));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'stopping_criteria'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'timestep_control'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'propagators'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'results_manager'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'sources'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'boundary_condition'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'convergency_study'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'visualizer'));
addpath(fullfile(SETTINGS.PROJECT_ROOT, 'console_logs'));

SETTINGS.SAVE_DIR = fullfile(SETTINGS.PROJECT_ROOT, 'figures_report_auto_v6');
SETTINGS.RUN_IMPLICIT = true;
SETTINGS.SAVE_FIGURES = true;
SETTINGS.SAVE_MAT = true;

SETTINGS.RUNTIME_REPEATS = 1;         % repetir tiempos para reducir ruido
SETTINGS.RUNTIME_DROP_FIRST = true;   % descarta primera corrida (warm-up)

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

SETTINGS.COST_REFERENCE_ORDER = 2.0; % integración temporal explícita con dt~h

if ~exist(SETTINGS.SAVE_DIR, 'dir')
    mkdir(SETTINGS.SAVE_DIR);
end

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
%  DATOS DE COSTE (ESTILO BENCHMARK)
% ==========================
n_cells_vec      = err_rus.ncells(:);
t_rusanov_repr   = err_rus.runtime(:);
t_hllc_repr      = err_hllc.runtime(:);
dt_rusanov_used  = cellfun(@(c) c.dt, rusanov_cases(:));
dt_hllc_used     = cellfun(@(c) c.dt, hllc_cases(:));
nsteps_rusanov   = ceil((SETTINGS.T_END - SETTINGS.T0) ./ dt_rusanov_used(:));
nsteps_hllc      = ceil((SETTINGS.T_END - SETTINGS.T0) ./ dt_hllc_used(:));

p_cost_rusanov = local_observed_orders(n_cells_vec, t_rusanov_repr);
p_cost_hllc    = local_observed_orders(n_cells_vec, t_hllc_repr);

coef_r = polyfit(log(n_cells_vec), log(t_rusanov_repr), 1);
coef_h = polyfit(log(n_cells_vec), log(t_hllc_repr), 1);

p_global_rusanov = coef_r(1);
p_global_hllc    = coef_h(1);

if SETTINGS.RUN_IMPLICIT && ~isempty(implicit_cases)
    imp_nodes_tag    = SETTINGS.IMPLICIT_MESH_LIST(:);
    imp_ncells       = cellfun(@(c) c.ncells,  implicit_cases(:));
    imp_h            = cellfun(@(c) c.h,       implicit_cases(:));
    imp_runtime      = cellfun(@(c) c.runtime, implicit_cases(:));
    imp_nnodes       = cellfun(@(c) c.nnodes,  implicit_cases(:));
    imp_nsteps       = ceil((SETTINGS.T_END - SETTINGS.T0) ./ cellfun(@(c) c.dt, implicit_cases(:)));
    imp_dt           = cellfun(@(c) c.dt,      implicit_cases(:));

    if numel(imp_ncells) >= 2
        p_cost_imp = local_observed_orders(imp_ncells(:), imp_runtime(:));
        coef_imp = polyfit(log(imp_ncells(:)), log(imp_runtime(:)), 1);
        p_global_imp = coef_imp(1);
    else
        p_cost_imp = [];
        coef_imp = [];
        p_global_imp = NaN;
    end
else
    imp_nodes_tag = [];
    imp_ncells    = [];
    imp_h         = [];
    imp_runtime   = [];
    imp_nnodes    = [];
    imp_nsteps    = [];
    imp_dt        = [];
    p_cost_imp    = [];
    coef_imp      = [];
    p_global_imp  = NaN;
end

%% =========================
%  FIGURA 4: COSTE COMPUTACIONAL (ESTILO BENCHMARK)
% ==========================
fprintf('\n=== TABLA RESUMEN DEL COSTE COMPUTACIONAL ===\n');

T_cost_exp = table( ...
    mesh_list(:), ...
    n_cells_vec(:), ...
    dt_rusanov_used(:), ...
    nsteps_rusanov(:), ...
    t_rusanov_repr(:), ...
    dt_hllc_used(:), ...
    nsteps_hllc(:), ...
    t_hllc_repr(:), ...
    t_hllc_repr(:) ./ t_rusanov_repr(:), ...
    [NaN; p_cost_rusanov(:)], ...
    [NaN; p_cost_hllc(:)], ...
    'VariableNames', { ...
        'nodes_tag', ...
        'num_cells', ...
        'dt_rusanov_s', ...
        'nsteps_rusanov', ...
        'time_rusanov_s', ...
        'dt_hllc_s', ...
        'nsteps_hllc', ...
        'time_hllc_s', ...
        'ratio_hllc_over_rusanov', ...
        'p_local_rusanov', ...
        'p_local_hllc' ...
    });

disp(T_cost_exp);

fprintf('\nOrden global Rusanov explícito = %.4f\n', p_global_rusanov);
fprintf('Orden global HLLC explícito    = %.4f\n', p_global_hllc);

fig = figure('Name', 'cost_explicit', 'Color', 'w');
loglog(n_cells_vec, t_rusanov_repr, '-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
loglog(n_cells_vec, t_hllc_repr,    '-s', 'LineWidth', 1.5, 'MarkerSize', 6);

grid on;
xlabel('Número de celdas');
ylabel('Tiempo de integración [s]');
title('Coste computacional de la integración temporal');

legend( ...
    sprintf('Rusanov exp. (p_{glob}=%.2f)', p_global_rusanov), ...
    sprintf('HLLC exp. (p_{glob}=%.2f)', p_global_hllc), ...
    'Location', 'northwest');

save_figure(fig, SETTINGS, 'fig04_coste_computacional');

%% =========================
%  FIGURA 5: ORDEN LOCAL DEL COSTE (ESTILO BENCHMARK)
% ==========================
fig = figure('Name', 'cost_order', 'Color', 'w');
plot(n_cells_vec(2:end), p_cost_rusanov, '-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
plot(n_cells_vec(2:end), p_cost_hllc,    '-s', 'LineWidth', 1.5, 'MarkerSize', 6);

yline(SETTINGS.COST_REFERENCE_ORDER, '--k', 'LineWidth', 1.2);
grid on;
xlabel('Número de celdas');
ylabel('Orden local observado');
title('Orden local del coste computacional de la integración temporal');

legend('Rusanov exp.', 'HLLC exp.', ...
    sprintf('Orden %.0f', SETTINGS.COST_REFERENCE_ORDER), ...
    'Location', 'best');

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
        'rich_rus', 'rich_hllc', 'gci_series_rus', 'gci_series_hllc', ...
        'T_cost_exp', ...
        'p_cost_rusanov', 'p_cost_hllc', 'p_cost_imp', ...
        'p_global_rusanov', 'p_global_hllc', 'p_global_imp');
end

fprintf('\nFiguras guardadas en:\n%s\n', SETTINGS.SAVE_DIR);

%% ========================================================================
%  FUNCIONES LOCALES
% ========================================================================

function out = simulate_case(nnodes, dt_value, flux_name, time_mode, SETTINGS)
    [nodes_file, cells_file, bc_files] = mesh_case_paths(nnodes);
    cells = mesh_processor(nodes_file, cells_file, bc_files);

    centroids_x = get_cell_centroid_x(cells);
    centroids_x = centroids_x(:);
    num_cells = numel(centroids_x);

    w = Config.INITIAL_CONDITIONS(centroids_x);
    t = SETTINGS.T0;

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
    for irun = 1:SETTINGS.RUNTIME_REPEATS
        w_run = w;
        t_run = SETTINGS.T0;

        tic;
        while t_run < SETTINGS.T_END - 1e-14
            dt = min(dt_value, SETTINGS.T_END - t_run);
            w_run = propagator(w_run, t_run, dt, problem);
            t_run = t_run + dt;
        end
        runtime_samples(irun) = toc;
    end

    if SETTINGS.RUNTIME_DROP_FIRST && numel(runtime_samples) >= 3
        runtime = median(runtime_samples(2:end));
    else
        runtime = median(runtime_samples);
    end

    w = w_run;
    t = t_run;

    rho  = w(1:num_cells);
    rhou = w(num_cells+1:2*num_cells);
    E    = w(2*num_cells+1:3*num_cells);
    u = rhou ./ rho;
    p = (SETTINGS.GAMMA - 1) .* (E - 0.5 .* (rhou.^2) ./ rho);
    ei = p ./ ((SETTINGS.GAMMA - 1) .* rho);

    [x_sorted, idx] = sort(centroids_x);
    rho = rho(idx);
    u = u(idx);
    p = p(idx);
    E = E(idx);

    out = struct();
    out.nnodes = nnodes;
    out.ncells = num_cells;
    out.h = mean(diff(x_sorted));
    out.dt = dt_value;
    out.t = t;
    out.runtime = runtime;
    out.runtime_samples = runtime_samples;
    out.x = x_sorted;
    out.rho = rho;
    out.u = u;
    out.p = p;
    out.E = E;
    out.ei = ei(idx);
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
        err.runtime(k) = cases_in{k}.runtime;
        err.rho(k) = mean(abs(cases_in{k}.rho - rho_ex));
        err.u(k)   = mean(abs(cases_in{k}.u - u_ex));
        err.p(k)   = mean(abs(cases_in{k}.p - p_ex));
    end
end

function ord = compute_observed_order(h, err)
    ord.h_mid = sqrt(h(1:end-1) .* h(2:end));
    ord.rho = abs(log(err.rho(2:end) ./ err.rho(1:end-1)) ./ log(h(2:end) ./ h(1:end-1)));
    ord.u   = abs(log(err.u(2:end) ./ err.u(1:end-1)) ./ log(h(2:end) ./ h(1:end-1)));
    ord.p   = abs(log(err.p(2:end) ./ err.p(1:end-1)) ./ log(h(2:end) ./ h(1:end-1)));
end

function p_local = local_observed_orders(N, t)
    N = N(:);
    t = t(:);

    if numel(N) ~= numel(t)
        error('N y t deben tener la misma longitud.');
    end
    if numel(N) < 2
        p_local = [];
        return;
    end

    p_local = zeros(length(N)-1,1);
    for i = 2:length(N)
        p_local(i-1) = log(t(i)/t(i-1)) / log(N(i)/N(i-1));
    end
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
        sprintf('%d nodos \\rightarrow malla gruesa', rich.medium.nnodes), ...
        sprintf('%d nodos \\rightarrow malla gruesa', rich.fine.nnodes), ...
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