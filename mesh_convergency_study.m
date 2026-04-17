clc; clear; close all;

addpath('convergency_study');

% =========================================================
% SETTINGS
% =========================================================
mesh_list = [66, 130, 258];
dt_value = 1e-5; % keep same dt for all meshes to reduce temporal contamination

% =========================================================
% RUN CASES
% =========================================================
cases = cell(numel(mesh_list), 1);

for k = 1:numel(mesh_list)
    fprintf('Running mesh %d nodes...\n', mesh_list(k));
    cases{k} = run_sod_case(mesh_list(k), dt_value);

    fprintf('  ncells   = %d\n', cases{k}.ncells);
    fprintf('  runtime  = %.6f s\n', cases{k}.runtime);
    fprintf('  t_final  = %.6f s\n', cases{k}.t);
end

coarse = cases{1}; % 66 nodes
medium = cases{2}; % 130 nodes
fine = cases{3}; % 258 nodes

% =========================================================
% PLOTS: SUPERPOSED SOLUTIONS
% =========================================================
figure;

subplot(3, 1, 1);
stairs(coarse.x, coarse.rho, 'LineWidth', 1.5); hold on;
stairs(medium.x, medium.rho, 'LineWidth', 1.5);
stairs(fine.x, fine.rho, 'LineWidth', 1.5);
grid on;
xlabel('x [m]');
ylabel('\rho [kg/m^3]');
title(sprintf('Density comparison at t \\approx %.4f s', fine.t));
legend('66 nodes', '130 nodes', '258 nodes', 'Location', 'best');

subplot(3, 1, 2);
stairs(coarse.x, coarse.u, 'LineWidth', 1.5); hold on;
stairs(medium.x, medium.u, 'LineWidth', 1.5);
stairs(fine.x, fine.u, 'LineWidth', 1.5);
grid on;
xlabel('x [m]');
ylabel('u [m/s]');
title('Velocity comparison');
legend('66 nodes', '130 nodes', '258 nodes', 'Location', 'best');

subplot(3, 1, 3);
stairs(coarse.x, coarse.p, 'LineWidth', 1.5); hold on;
stairs(medium.x, medium.p, 'LineWidth', 1.5);
stairs(fine.x, fine.p, 'LineWidth', 1.5);
grid on;
xlabel('x [m]');
ylabel('p [Pa]');
title('Pressure comparison');
legend('66 nodes', '130 nodes', '258 nodes', 'Location', 'best');

% =========================================================
% COST ANALYSIS
% =========================================================
ncells = [coarse.ncells, medium.ncells, fine.ncells];
runtimes = [coarse.runtime, medium.runtime, fine.runtime];
nsteps = [coarse.nsteps_est, medium.nsteps_est, fine.nsteps_est];
work_units = ncells .* nsteps;

figure;
plot(ncells, runtimes, 'o-', 'LineWidth', 2);
grid on;
xlabel('Number of cells');
ylabel('Runtime [s]');
title('Computational cost vs mesh size');

cost_table = table( ...
    mesh_list(:), ...
    ncells(:), ...
    nsteps(:), ...
    runtimes(:), ...
    work_units(:), ...
    'VariableNames', {'NodesFileTag', 'Ncells', 'EstimatedSteps', 'Runtime_s', 'CellUpdates'} ...
);

disp('=== COST TABLE ===');
disp(cost_table);

% =========================================================
% RICHARDSON + GCI
% =========================================================
fprintf('\n=== RICHARDSON / GCI ANALYSIS ===\n');

rho_metrics = richardson_gci_three_grids(coarse.x, coarse.rho, medium.x, medium.rho, fine.x, fine.rho);
u_metrics = richardson_gci_three_grids(coarse.x, coarse.u, medium.x, medium.u, fine.x, fine.u);
p_metrics = richardson_gci_three_grids(coarse.x, coarse.p, medium.x, medium.p, fine.x, fine.p);

metrics_table = table( ...
    ["density"; "velocity"; "pressure"], ...
    [rho_metrics.p; u_metrics.p; p_metrics.p], ...
    [rho_metrics.E32; u_metrics.E32; p_metrics.E32], ...
    [rho_metrics.E21; u_metrics.E21; p_metrics.E21], ...
    [rho_metrics.GCI32; u_metrics.GCI32; p_metrics.GCI32], ...
    [rho_metrics.GCI21; u_metrics.GCI21; p_metrics.GCI21], ...
    [rho_metrics.asymptotic_ratio; u_metrics.asymptotic_ratio; p_metrics.asymptotic_ratio], ...
    'VariableNames', {'Variable', 'ObservedOrder_p', 'E32', 'E21', 'GCI32_percent', 'GCI21_percent', 'AsymptoticRatio'} ...
);

disp(metrics_table);

% Optional: extrapolated fine solution for density on coarse grid
figure;
stairs(coarse.x, rho_metrics.phi3, 'LineWidth', 1.5); hold on;
stairs(coarse.x, rho_metrics.phi2_on_3, 'LineWidth', 1.5);
stairs(coarse.x, rho_metrics.phi1_on_3, 'LineWidth', 1.5);
plot(coarse.x, rho_metrics.phi_ext, '--', 'LineWidth', 2);
grid on;
xlabel('x [m]');
ylabel('\rho [kg/m^3]');
title('Density: coarse / medium / fine / Richardson extrapolated');
legend('coarse', 'medium->coarse', 'fine->coarse', 'Richardson extrapolated', 'Location', 'best');
