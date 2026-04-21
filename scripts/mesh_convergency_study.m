clc; clear; close all;

addpath('convergency_study');
addpath('propagators');

%% =========================================================
% SETTINGS
% =========================================================
mesh_list = [66, 130, 258, 514, 1026, 2050];

% Queremos dt ~ h para que el coste computacional tienda a orden 2
% Tomamos como referencia el dt de la malla más fina
dt_finest = 1e-5;

prop_explicit = @(w, t, dt, f) fw_euler_sparse(w, t, dt, f);
prop_implicit = @(w, t, dt, f) bw_euler_sparse(w, t, dt, f);

nMeshes = numel(mesh_list);

% Estimación previa de celdas a partir del tag de malla
% (66->32, 130->64, ..., 4098->2048)
ncells_nominal = (mesh_list - 2) / 2;
h_nominal = 1 ./ ncells_nominal;

% dt escalado con h para mantener algo equivalente a CFL constante
dt_list = dt_finest * (h_nominal / h_nominal(end));

fprintf('=== DT POR MALLA (ESCALADO CON h) ===\n');
for k = 1:nMeshes
    fprintf('Mesh %d nodes -> ncells ~ %d, h ~ %.4e, dt = %.4e\n', ...
        mesh_list(k), ncells_nominal(k), h_nominal(k), dt_list(k));
end

%% =========================================================
% RUN CASOS - EXPLÍCITO
% =========================================================
fprintf('\n=== MÉTODO EXPLÍCITO ===\n');
cases_expl = cell(nMeshes, 1);

for k = 1:nMeshes
    fprintf('Running mesh %d nodes...\n', mesh_list(k));
    cases_expl{k} = run_sod_case(mesh_list(k), dt_list(k), prop_explicit);
    fprintf('  ncells  = %d\n', cases_expl{k}.ncells);
    fprintf('  runtime = %.6f s\n', cases_expl{k}.runtime);
    fprintf('  dt      = %.4e\n', dt_list(k));
    fprintf('  t_final = %.6f s\n', cases_expl{k}.t);
end

%% =========================================================
% RUN CASOS - IMPLÍCITO
% =========================================================
fprintf('\n=== MÉTODO IMPLÍCITO ===\n');
cases_impl = cell(nMeshes, 1);

for k = 1:nMeshes
    fprintf('Running mesh %d nodes...\n', mesh_list(k));
    cases_impl{k} = run_sod_case(mesh_list(k), dt_list(k), prop_implicit);
    fprintf('  ncells  = %d\n', cases_impl{k}.ncells);
    fprintf('  runtime = %.6f s\n', cases_impl{k}.runtime);
    fprintf('  dt      = %.4e\n', dt_list(k));
    fprintf('  t_final = %.6f s\n', cases_impl{k}.t);
end

%% =========================================================
% EXTRAER DATOS
% =========================================================
ncells = cellfun(@(c) c.ncells, cases_expl);
h_vals = 1 ./ ncells;

runtimes_expl = cellfun(@(c) c.runtime, cases_expl);
runtimes_impl = cellfun(@(c) c.runtime, cases_impl);

nsteps_expl = cellfun(@(c) c.nsteps_est, cases_expl);
nsteps_impl = cellfun(@(c) c.nsteps_est, cases_impl);

work_units_expl = ncells .* nsteps_expl;
work_units_impl = ncells .* nsteps_impl;

% Para Richardson/GCI usamos las 3 mallas más finas del método explícito
coarse = cases_expl{end-2};
medium = cases_expl{end-1};
fine   = cases_expl{end};

%% =========================================================
% ÍNDICES PARA FIGURA 1
% Representar todas excepto 1026 y 2050, pero sí incluir 4098
%% =========================================================
plot_idx = ~ismember(mesh_list, [1026, 2050]);

%% =========================================================
% FIGURA 1: SOLUCIONES SUPERPUESTAS
%% =========================================================
figure;

subplot(3,1,1); hold on;
for k = find(plot_idx)
    stairs(cases_expl{k}.x, cases_expl{k}.rho, 'LineWidth', 1.5);
end
grid on;
xlabel('x [m]');
ylabel('\rho [kg/m^3]');
title(sprintf('Density comparison at t \\approx %.4f s', cases_expl{end}.t));
legend(compose('%d nodes', mesh_list(plot_idx)), 'Location', 'best');

subplot(3,1,2); hold on;
for k = find(plot_idx)
    stairs(cases_expl{k}.x, cases_expl{k}.u, 'LineWidth', 1.5);
end
grid on;
xlabel('x [m]');
ylabel('u [m/s]');
title('Velocity comparison');
legend(compose('%d nodes', mesh_list(plot_idx)), 'Location', 'best');

subplot(3,1,3); hold on;
for k = find(plot_idx)
    stairs(cases_expl{k}.x, cases_expl{k}.p, 'LineWidth', 1.5);
end
grid on;
xlabel('x [m]');
ylabel('p [Pa]');
title('Pressure comparison');
legend(compose('%d nodes', mesh_list(plot_idx)), 'Location', 'best');

%% =========================================================
% FIGURA 2: COSTE COMPUTACIONAL
%% =========================================================
figure;
plot(ncells, runtimes_expl, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
plot(ncells, runtimes_impl, 'ro-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
xlabel('Number of cells');
ylabel('Runtime [s]');
title('Computational cost vs mesh size');
legend('Explicit', 'Implicit', 'Location', 'best');

%% =========================================================
% TABLA DE COSTES
%% =========================================================
cost_table = table( ...
    mesh_list(:), ...
    ncells(:), ...
    h_vals(:), ...
    dt_list(:), ...
    nsteps_expl(:), ...
    nsteps_impl(:), ...
    runtimes_expl(:), ...
    runtimes_impl(:), ...
    work_units_expl(:), ...
    work_units_impl(:), ...
    'VariableNames', { ...
    'NodesFileTag', ...
    'Ncells', ...
    'h_m', ...
    'dt_used_s', ...
    'EstimatedSteps_Expl', ...
    'EstimatedSteps_Impl', ...
    'Runtime_Expl_s', ...
    'Runtime_Impl_s', ...
    'CellUpdates_Expl', ...
    'CellUpdates_Impl'} ...
);

fprintf('\n=== COST TABLE ===\n');
disp(cost_table);

%% =========================================================
% FIGURA 3: ORDEN DEL COSTE COMPUTACIONAL
%% =========================================================
order_expl = abs(log(runtimes_expl(2:end) ./ runtimes_expl(1:end-1)) ./ ...
                 log(h_vals(2:end) ./ h_vals(1:end-1)));

order_impl = abs(log(runtimes_impl(2:end) ./ runtimes_impl(1:end-1)) ./ ...
                 log(h_vals(2:end) ./ h_vals(1:end-1)));

h_mid = sqrt(h_vals(1:end-1) .* h_vals(2:end));

fprintf('\n=== ORDEN DEL COSTE COMPUTACIONAL ===\n');
for k = 1:numel(h_mid)
    fprintf('h = %.2e | explícito: %.2f | implícito: %.2f\n', ...
        h_mid(k), order_expl(k), order_impl(k));
end

h_plot = flip(h_mid);
order_expl_plot = flip(order_expl);
order_impl_plot = flip(order_impl);

figure;
semilogx(h_plot, order_expl_plot, 'bo-', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
hold on;
semilogx(h_plot, order_impl_plot, 'ro-', 'LineWidth', 2, ...
    'MarkerFaceColor', 'r', 'MarkerSize', 8);
yline(2, 'k--', 'LineWidth', 1.5);
grid on;
xlabel('h [m]');
ylabel('Orden del método [-]');
title('Evolución del orden del coste computacional');
legend('Método explícito', 'Método implícito', 'Orden 2 (referencia)', ...
    'Location', 'best');

y_max = max([order_expl(:); order_impl(:); 2.2]);
ylim([0, y_max * 1.1]);

%% =========================================================
% TABLA ORDEN DEL COSTE COMPUTACIONAL
%% =========================================================
order_table = table( ...
    h_mid(:), ...
    order_expl(:), ...
    order_impl(:), ...
    'VariableNames', {'h_mid_m', 'Order_Explicit', 'Order_Implicit'} ...
);

fprintf('\n=== ORDER TABLE ===\n');
disp(order_table);

%% =========================================================
% RICHARDSON + GCI (3 MALLAS MÁS FINAS, MÉTODO EXPLÍCITO)
%% =========================================================
fprintf('\n=== RICHARDSON / GCI ANALYSIS ===\n');
fprintf('Using meshes for Richardson/GCI: %d, %d, %d nodes\n', ...
    mesh_list(end-2), mesh_list(end-1), mesh_list(end));

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

%% =========================================================
% FIGURA 4: RICHARDSON EXTRAPOLATION - DENSITY
%% =========================================================
figure;
stairs(coarse.x, rho_metrics.phi3, 'LineWidth', 1.5); hold on;
stairs(coarse.x, rho_metrics.phi2_on_3, 'LineWidth', 1.5);
stairs(coarse.x, rho_metrics.phi1_on_3, 'LineWidth', 1.5);
plot(coarse.x, rho_metrics.phi_ext, '--', 'LineWidth', 2);
grid on;
xlabel('x [m]');
ylabel('\rho [kg/m^3]');
title(sprintf('Density: %d / %d / %d nodes + Richardson extrapolated', ...
    mesh_list(end-2), mesh_list(end-1), mesh_list(end)));
legend( ...
    sprintf('%d nodes', mesh_list(end-2)), ...
    sprintf('%d nodes -> coarse grid', mesh_list(end-1)), ...
    sprintf('%d nodes -> coarse grid', mesh_list(end)), ...
    'Richardson extrapolated', ...
    'Location', 'best');

%% =========================================================
% FIGURA 5: CONVERGENCIA DE MALLA (ENERGÍA CINÉTICA ESPECÍFICA MEDIA)
%% =========================================================
Ek_mean = zeros(nMeshes,1);

for k = 1:nMeshes
    u_k = cases_expl{k}.u(:);
    Ek_local = 0.5 * (u_k.^2);
    Ek_mean(k) = mean(Ek_local);
end

figure;
semilogx(h_vals, Ek_mean, 'o-', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
grid on;
xlabel('h (m)');
ylabel('Energía cinética (J/kg)');
title('Convergencia de malla');

%% =========================================================
% TABLA DE CONVERGENCIA DE MALLA
%% =========================================================
mesh_conv_table = table( ...
    mesh_list(:), ...
    ncells(:), ...
    h_vals(:), ...
    Ek_mean(:), ...
    'VariableNames', {'NodesFileTag', 'Ncells', 'h_m', 'MeanKineticEnergy'} ...
);

fprintf('\n=== MESH CONVERGENCE TABLE ===\n');
disp(mesh_conv_table);

%% =========================================================
% GCI / AR POR TRÍOS CONSECUTIVOS USANDO h Y Ek_mean
%% =========================================================
fprintf('\n=== GCI / AR POR TRÍOS CONSECUTIVOS ===\n');

h = h_vals(:);
Ek = Ek_mean(:);
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

resultados_gci_ar = table('Size',[nTrios 10], ...
    'VariableTypes', {'double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'h3','h2','h1','r32','r21','s','p','GCI32','GCI21','AR'});

for i = 1:nTrios
    h3 = h(i);
    h2 = h(i+1);
    h1 = h(i+2);

    Ek3 = Ek(i);
    Ek2 = Ek(i+1);
    Ek1 = Ek(i+2);

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

    ea32 = abs((Ek2 - Ek3) / Ek2);
    ea21 = abs((Ek1 - Ek2) / Ek1);

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

    resultados_gci_ar{i,:} = [h3, h2, h1, r32, r21, s, p, GCI32, GCI21, AR];
end

disp('=========================================================')
disp('RESULTADOS GCI / AR')
disp('=========================================================')
disp(resultados_gci_ar)

%% =========================================================
% FIGURA 6: GRÁFICA DEL AR
%% =========================================================
figure;
plot(1:nTrios, resultados_gci_ar.AR, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
hold on;
yline(1, '--', 'LineWidth', 1.2);
grid on;
xlabel('Índice del trío de mallas');
ylabel('AR');
title('Rango asintótico (AR) por trío de mallas');

labels = strings(nTrios,1);
for i = 1:nTrios
    labels(i) = sprintf('[%.4g, %.4g, %.4g]', ...
        resultados_gci_ar.h3(i), resultados_gci_ar.h2(i), resultados_gci_ar.h1(i));
end
xticks(1:nTrios);
xticklabels(labels);
xtickangle(20);

%% =========================================================
% FIGURA 7: GRÁFICA DEL GCI
%% =========================================================
figure;
plot(1:nTrios, 100*resultados_gci_ar.GCI32, 's-', 'LineWidth', 1.8, 'MarkerSize', 7);
hold on;
plot(1:nTrios, 100*resultados_gci_ar.GCI21, 'o-', 'LineWidth', 1.8, 'MarkerSize', 7);
grid on;
xlabel('Índice del trío de mallas');
ylabel('GCI [%]');
title('Índice de convergencia de malla');
legend('GCI_{32}', 'GCI_{21}', 'Location', 'best');

%% =========================================================
% MENSAJE FINAL
%% =========================================================
disp(' ')
disp('Interpretación:')
disp('- Si AR ≈ 1, el trío de mallas está en régimen asintótico.')
disp('- Si AR se aleja claramente de 1, conviene refinar más la malla.')
