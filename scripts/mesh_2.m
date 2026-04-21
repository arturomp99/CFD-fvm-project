clc; clear; close all;

addpath('convergency_study');
addpath('propagators');
addpath('mesh_processing');
addpath(genpath('utils'));
addpath('constants');
addpath('initial_conditions');
addpath('problems/fvm_1D_euler');
addpath(genpath('convective_flux'));
addpath('stopping_criteria');
addpath('timestep_control');
addpath('results_manager');
addpath('boundary_condition');
addpath('sources');

%% =========================================================
% SETTINGS
% =========================================================
mesh_list = [66, 130, 258, 514, 1026, 2050, 4098];

% Queremos dt ~ h
dt_finest = 1e-5;

prop_explicit = @(w, t, dt, f) fw_euler_sparse(w, t, dt, f);

nMeshes = numel(mesh_list);

% Métodos explícitos a comparar
method_names = { ...
                    'Rusanov', ...
                    'HLLC' ...
                };

flux_methods = { ...
                    @(state, cells) rusanov_interpolator(state, cells), ...
                    @(state, cells) hllc_interpolator(state, cells) ...
                };

nMethods = numel(flux_methods);

% Estimación previa de celdas a partir del tag de malla
ncells_nominal = (mesh_list - 2) / 2;
h_nominal = 1 ./ ncells_nominal;

% dt escalado con h
dt_list = dt_finest * (h_nominal / h_nominal(end));

fprintf('=== DT POR MALLA (ESCALADO CON h) ===\n');

for k = 1:nMeshes
    fprintf('Mesh %d nodes -> ncells ~ %d, h ~ %.4e, dt = %.4e\n', ...
        mesh_list(k), ncells_nominal(k), h_nominal(k), dt_list(k));
end

%% =========================================================
% RUN CASOS - MÉTODOS EXPLÍCITOS
% =========================================================
cases_explicit = cell(nMethods, nMeshes);

for m = 1:nMethods
    fprintf('\n=====================================================\n');
    fprintf('MÉTODO EXPLÍCITO: %s\n', method_names{m});
    fprintf('=====================================================\n');

    for k = 1:nMeshes
        fprintf('Running mesh %d nodes...\n', mesh_list(k));

        cases_explicit{m, k} = run_sod_case_custom_flux( ...
            mesh_list(k), ...
            dt_list(k), ...
            prop_explicit, ...
            flux_methods{m} ...
        );

        fprintf('  ncells  = %d\n', cases_explicit{m, k}.ncells);
        fprintf('  runtime = %.6f s\n', cases_explicit{m, k}.runtime);
        fprintf('  dt      = %.4e\n', dt_list(k));
        fprintf('  t_final = %.6f s\n', cases_explicit{m, k}.t);
    end

end

%% =========================================================
% EXTRAER DATOS
% =========================================================
ncells = cellfun(@(c) c.ncells, cases_explicit(1, :));
h_vals = 1 ./ ncells;

runtimes = zeros(nMethods, nMeshes);
nsteps = zeros(nMethods, nMeshes);
work_units = zeros(nMethods, nMeshes);

for m = 1:nMethods
    runtimes(m, :) = cellfun(@(c) c.runtime, cases_explicit(m, :));
    nsteps(m, :) = cellfun(@(c) c.nsteps_est, cases_explicit(m, :));
    work_units(m, :) = ncells .* nsteps(m, :);
end

%% =========================================================
% ÍNDICES PARA FIGURA 1
% Representar todas excepto 1026 y 2050, pero sí incluir 4098
%% =========================================================
plot_idx = ~ismember(mesh_list, [1026, 2050]);

%% =========================================================
% FIGURA 1: SOLUCIONES SUPERPUESTAS PARA LA MALLA MÁS FINA
%% =========================================================
kfine = nMeshes;

figure;

subplot(3, 1, 1); hold on;

for m = 1:nMethods
    stairs(cases_explicit{m, kfine}.x, cases_explicit{m, kfine}.rho, 'LineWidth', 1.5);
end

grid on;
xlabel('x [m]');
ylabel('\rho [kg/m^3]');
title(sprintf('Density comparison at t \\approx %.4f s (%d nodes)', ...
    cases_explicit{1, kfine}.t, mesh_list(kfine)));
legend(method_names, 'Location', 'best');

subplot(3, 1, 2); hold on;

for m = 1:nMethods
    stairs(cases_explicit{m, kfine}.x, cases_explicit{m, kfine}.u, 'LineWidth', 1.5);
end

grid on;
xlabel('x [m]');
ylabel('u [m/s]');
title('Velocity comparison');
legend(method_names, 'Location', 'best');

subplot(3, 1, 3); hold on;

for m = 1:nMethods
    stairs(cases_explicit{m, kfine}.x, cases_explicit{m, kfine}.p, 'LineWidth', 1.5);
end

grid on;
xlabel('x [m]');
ylabel('p [Pa]');
title('Pressure comparison');
legend(method_names, 'Location', 'best');

%% =========================================================
% FIGURA 2: SOLUCIONES SUPERPUESTAS POR MALLA PARA HLLC
%% =========================================================
% Escogemos HLLC como método de referencia para convergencia espacial
ref_method_idx = find(strcmp(method_names, 'HLLC'), 1);

if isempty(ref_method_idx)
    ref_method_idx = 1;
end

figure;

subplot(3, 1, 1); hold on;

for k = find(plot_idx)
    stairs(cases_explicit{ref_method_idx, k}.x, cases_explicit{ref_method_idx, k}.rho, 'LineWidth', 1.5);
end

grid on;
xlabel('x [m]');
ylabel('\rho [kg/m^3]');
title(sprintf('Density comparison (%s)', method_names{ref_method_idx}));
legend(compose('%d nodes', mesh_list(plot_idx)), 'Location', 'best');

subplot(3, 1, 2); hold on;

for k = find(plot_idx)
    stairs(cases_explicit{ref_method_idx, k}.x, cases_explicit{ref_method_idx, k}.u, 'LineWidth', 1.5);
end

grid on;
xlabel('x [m]');
ylabel('u [m/s]');
title('Velocity comparison');
legend(compose('%d nodes', mesh_list(plot_idx)), 'Location', 'best');

subplot(3, 1, 3); hold on;

for k = find(plot_idx)
    stairs(cases_explicit{ref_method_idx, k}.x, cases_explicit{ref_method_idx, k}.p, 'LineWidth', 1.5);
end

grid on;
xlabel('x [m]');
ylabel('p [Pa]');
title('Pressure comparison');
legend(compose('%d nodes', mesh_list(plot_idx)), 'Location', 'best');

%% =========================================================
% FIGURA 3: COSTE COMPUTACIONAL
%% =========================================================
figure; hold on;

for m = 1:nMethods
    plot(ncells, runtimes(m, :), 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
end

grid on;
xlabel('Number of cells');
ylabel('Runtime [s]');
title('Computational cost vs mesh size');
legend(method_names, 'Location', 'best');

%% =========================================================
% TABLA DE COSTES
%% =========================================================
for m = 1:nMethods
    cost_table = table( ...
        repmat(string(method_names{m}), nMeshes, 1), ...
        mesh_list(:), ...
        ncells(:), ...
        h_vals(:), ...
        dt_list(:), ...
        nsteps(m, :).', ...
        runtimes(m, :).', ...
        work_units(m, :).', ...
        'VariableNames', { ...
           'Method', ...
           'NodesFileTag', ...
           'Ncells', ...
           'h_m', ...
           'dt_used_s', ...
           'EstimatedSteps', ...
           'Runtime_s', ...
       'CellUpdates'} ...
    );

    fprintf('\n=== COST TABLE: %s ===\n', method_names{m});
    disp(cost_table);
end

%% =========================================================
% FIGURA 4: ORDEN DEL COSTE COMPUTACIONAL
%% =========================================================
figure; hold on;

for m = 1:nMethods
    order_m = abs(log(runtimes(m, 2:end) ./ runtimes(m, 1:end - 1)) ./ ...
        log(h_vals(2:end) ./ h_vals(1:end - 1)));

    h_mid = sqrt(h_vals(1:end - 1) .* h_vals(2:end));
    h_plot = flip(h_mid);
    order_plot = flip(order_m);

    semilogx(h_plot, order_plot, 'o-', 'LineWidth', 2, ...
        'MarkerFaceColor', 'auto', 'MarkerSize', 8);
end

yline(2, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('h [m]');
ylabel('Orden del coste [-]');
title('Evolución del orden del coste computacional');
legend([method_names, {'Orden 2 (referencia)'}], 'Location', 'best');

%% =========================================================
% TABLA ORDEN DEL COSTE COMPUTACIONAL
%% =========================================================
for m = 1:nMethods
    order_m = abs(log(runtimes(m, 2:end) ./ runtimes(m, 1:end - 1)) ./ ...
        log(h_vals(2:end) ./ h_vals(1:end - 1)));
    h_mid = sqrt(h_vals(1:end - 1) .* h_vals(2:end));

    order_table = table( ...
        h_mid(:), ...
        order_m(:), ...
        'VariableNames', {'h_mid_m', 'Order'} ...
    );

    fprintf('\n=== ORDER TABLE: %s ===\n', method_names{m});
    disp(order_table);
end

%% =========================================================
% RICHARDSON + GCI (3 MALLAS MÁS FINAS, PARA CADA MÉTODO)
%% =========================================================
fprintf('\n=== RICHARDSON / GCI ANALYSIS ===\n');

for m = 1:nMethods
    fprintf('\n---------------------------------------------\n');
    fprintf('Método: %s\n', method_names{m});
    fprintf('Using meshes: %d, %d, %d nodes\n', ...
        mesh_list(end - 2), mesh_list(end - 1), mesh_list(end));

    coarse = cases_explicit{m, end - 2};
    medium = cases_explicit{m, end - 1};
    fine = cases_explicit{m, end};

    rho_metrics = richardson_gci_three_grids( ...
        coarse.x, coarse.rho, medium.x, medium.rho, fine.x, fine.rho);

    u_metrics = richardson_gci_three_grids( ...
        coarse.x, coarse.u, medium.x, medium.u, fine.x, fine.u);

    p_metrics = richardson_gci_three_grids( ...
        coarse.x, coarse.p, medium.x, medium.p, fine.x, fine.p);

    metrics_table = table( ...
        ["density"; "velocity"; "pressure"], ...
        [rho_metrics.p; u_metrics.p; p_metrics.p], ...
        [rho_metrics.E32; u_metrics.E32; p_metrics.E32], ...
        [rho_metrics.E21; u_metrics.E21; p_metrics.E21], ...
        [rho_metrics.GCI32; u_metrics.GCI32; p_metrics.GCI32], ...
        [rho_metrics.GCI21; u_metrics.GCI21; p_metrics.GCI21], ...
        [rho_metrics.asymptotic_ratio; u_metrics.asymptotic_ratio; p_metrics.asymptotic_ratio], ...
        'VariableNames', { ...
           'Variable', ...
           'ObservedOrder_p', ...
           'E32', ...
           'E21', ...
           'GCI32_percent', ...
           'GCI21_percent', ...
       'AsymptoticRatio'} ...
    );

    disp(metrics_table);

    % Figura de densidad extrapolada
    figure;
    stairs(coarse.x, rho_metrics.phi3, 'LineWidth', 1.5); hold on;
    stairs(coarse.x, rho_metrics.phi2_on_3, 'LineWidth', 1.5);
    stairs(coarse.x, rho_metrics.phi1_on_3, 'LineWidth', 1.5);
    plot(coarse.x, rho_metrics.phi_ext, '--', 'LineWidth', 2);
    grid on;
    xlabel('x [m]');
    ylabel('\rho [kg/m^3]');
    title(sprintf('Density - %s: %d / %d / %d nodes + Richardson extrapolated', ...
        method_names{m}, mesh_list(end - 2), mesh_list(end - 1), mesh_list(end)));
    legend( ...
        sprintf('%d nodes', mesh_list(end - 2)), ...
        sprintf('%d nodes -> coarse grid', mesh_list(end - 1)), ...
        sprintf('%d nodes -> coarse grid', mesh_list(end)), ...
        'Richardson extrapolated', ...
        'Location', 'best');
end

%% =========================================================
% FIGURA 5: CONVERGENCIA DE MALLA (ENERGÍA CINÉTICA ESPECÍFICA MEDIA)
%% =========================================================
Ek_mean = zeros(nMethods, nMeshes);

for m = 1:nMethods

    for k = 1:nMeshes
        u_k = cases_explicit{m, k}.u(:);
        Ek_local = 0.5 * (u_k .^ 2);
        Ek_mean(m, k) = mean(Ek_local);
    end

end

figure; hold on;

for m = 1:nMethods
    semilogx(h_vals, Ek_mean(m, :), 'o-', 'LineWidth', 2, ...
        'MarkerFaceColor', 'auto', 'MarkerSize', 8);
end

grid on;
xlabel('h [m]');
ylabel('Energía cinética específica media [J/kg]');
title('Convergencia de malla');
legend(method_names, 'Location', 'best');

%% =========================================================
% TABLAS DE CONVERGENCIA DE MALLA
%% =========================================================
for m = 1:nMethods
    mesh_conv_table = table( ...
        repmat(string(method_names{m}), nMeshes, 1), ...
        mesh_list(:), ...
        ncells(:), ...
        h_vals(:), ...
        Ek_mean(m, :).', ...
        'VariableNames', {'Method', 'NodesFileTag', 'Ncells', 'h_m', 'MeanKineticEnergy'} ...
    );

    fprintf('\n=== MESH CONVERGENCE TABLE: %s ===\n', method_names{m});
    disp(mesh_conv_table);
end

%% =========================================================
% GCI / AR POR TRÍOS CONSECUTIVOS USANDO h Y Ek_mean
%% =========================================================
fprintf('\n=== GCI / AR POR TRÍOS CONSECUTIVOS ===\n');

for m = 1:nMethods
    fprintf('\n---------------------------------------------\n');
    fprintf('Método: %s\n', method_names{m});

    h = h_vals(:);
    Ek = Ek_mean(m, :).';
    Fs = 1.25;

    if numel(h) ~= numel(Ek)
        error('h y Ek deben tener la misma longitud.');
    end

    if numel(h) < 3
        error('Necesitas al menos 3 mallas para calcular GCI y AR.');
    end

    if any(diff(h) >= 0)
        error('Los valores de h deben ir en orden estrictamente decreciente: h(1)>h(2)>h(3)>...');
    end

    n = numel(h);
    nTrios = n - 2;

    resultados_gci_ar = table('Size', [nTrios 10], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'h3', 'h2', 'h1', 'r32', 'r21', 's', 'p', 'GCI32', 'GCI21', 'AR'});

    for i = 1:nTrios
        h3 = h(i);
        h2 = h(i + 1);
        h1 = h(i + 2);

        Ek3 = Ek(i);
        Ek2 = Ek(i + 1);
        Ek1 = Ek(i + 2);

        r32 = h3 / h2;
        r21 = h2 / h1;

        e32 = Ek3 - Ek2;
        e21 = Ek2 - Ek1;

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
            p0 = abs(log(abs(e32 / e21)) / log(r21));

            if ~isfinite(p0) || p0 == 0
                p0 = 1;
            end

            p = p0;

            for kiter = 1:100
                num = r21 ^ p - s;
                den = r32 ^ p - s;

                if abs(den) < eps || abs(num) < eps
                    p = NaN;
                    break;
                end

                p_new = abs((log(abs(e32 / e21)) + log(abs(num / den))) / log(r21));

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

        ea32 = abs((Ek2 - Ek3) / Ek2);
        ea21 = abs((Ek1 - Ek2) / Ek1);

        if isnan(p) || abs(r32 ^ p - 1) < eps || abs(r21 ^ p - 1) < eps
            GCI32 = NaN;
            GCI21 = NaN;
            AR = NaN;
        else
            GCI32 = Fs * ea32 / (r32 ^ p - 1);
            GCI21 = Fs * ea21 / (r21 ^ p - 1);

            if abs(GCI32) < eps
                AR = NaN;
            else
                AR = (GCI21 * r21 ^ p) / GCI32;
            end

        end

        resultados_gci_ar{i, :} = [h3, h2, h1, r32, r21, s, p, GCI32, GCI21, AR];
    end

    disp(resultados_gci_ar)

    figure;
    plot(1:nTrios, resultados_gci_ar.AR, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
    hold on;
    yline(1, '--', 'LineWidth', 1.2);
    grid on;
    xlabel('Índice del trío de mallas');
    ylabel('AR');
    title(sprintf('Rango asintótico (AR) por trío de mallas - %s', method_names{m}));

    labels = strings(nTrios, 1);

    for i = 1:nTrios
        labels(i) = sprintf('[%.4g, %.4g, %.4g]', ...
            resultados_gci_ar.h3(i), resultados_gci_ar.h2(i), resultados_gci_ar.h1(i));
    end

    xticks(1:nTrios);
    xticklabels(labels);
    xtickangle(20);

    figure;
    plot(1:nTrios, 100 * resultados_gci_ar.GCI32, 's-', 'LineWidth', 1.8, 'MarkerSize', 7);
    hold on;
    plot(1:nTrios, 100 * resultados_gci_ar.GCI21, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
    grid on;
    xlabel('Índice del trío de mallas');
    ylabel('GCI [%]');
    title(sprintf('Índice de convergencia de malla - %s', method_names{m}));
    legend('GCI_{32}', 'GCI_{21}', 'Location', 'best');
end

%% =========================================================
% MENSAJE FINAL
%% =========================================================
disp(' ')
disp('Interpretación:')
disp('- Aquí comparas distintos métodos de flujo explícitos con el mismo integrador temporal.')
disp('- Si AR ≈ 1, el trío de mallas está en régimen asintótico.')
disp('- Si AR se aleja claramente de 1, conviene refinar más la malla.')

%% =========================================================
% FUNCIONES LOCALES
%% =========================================================
function out = run_sod_case_custom_flux(nnodes, dt_value, propagator_in, flux_handle)

    [nodes_file, cells_file, bc_files] = mesh_case_paths(nnodes);

    cells = mesh_processor(nodes_file, cells_file, bc_files);
    num_cells = length(cells);

    centroids_x = get_cell_centroid_x(cells);
    centroids_x = centroids_x(:);

    w0 = Config.INITIAL_CONDITIONS(centroids_x);

    problem = @(state, time) fvm_1D_euler_custom(state, cells, time, flux_handle);

    timestep_calculator = @(w, t) constant_dt(w, t, dt_value);

    stopping_condition = Config.STOPPING_CONDITION;
    manager = Config.RESULTS_MANAGER;

    tic;
    results = solver( ...
        w0, ...
        Config.T0, ...
        problem, ...
        propagator_in, ...
        timestep_calculator, ...
        stopping_condition, ...
        manager ...
    );
    runtime = toc;

    if isempty(results)
        error('run_sod_case_custom_flux:EmptyResults', ...
        'The solver returned an empty results array.');
    end

    t_vec = results(:, 1);
    [~, idx_t] = min(abs(t_vec - Config.T_END));
    t_plot = t_vec(idx_t);

    rho = results(idx_t, 2:num_cells + 1).';
    rhou = results(idx_t, num_cells + 2:2 * num_cells + 1).';
    E = results(idx_t, 2 * num_cells + 2:3 * num_cells + 1).';

    u = rhou ./ rho;
    p = (Air.GAMMA - 1) .* (E - 0.5 .* (rhou .^ 2) ./ rho);

    [x_sorted, idx_sort] = sort(centroids_x);

    out = struct();
    out.nnodes = nnodes;
    out.ncells = num_cells;
    out.x = x_sorted;
    out.rho = rho(idx_sort);
    out.u = u(idx_sort);
    out.p = p(idx_sort);
    out.t = t_plot;
    out.runtime = runtime;
    out.dt = dt_value;
    out.nsteps_est = round((Config.T_END - Config.T0) / dt_value);
    out.results = results;
end

function [A, b] = fvm_1D_euler_custom(state, cells, t, flux_handle)
    [A, b] = flux_handle(state, cells);

    centroids_x = get_cell_centroid_x(cells);
    [A_source, b_source] = Config.SOURCE_TERMS(state, centroids_x, t);

    A = A + A_source;
    b = b + b_source;
end
