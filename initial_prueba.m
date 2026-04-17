clc; clear; close all;

addpath('mesh_processing');
addpath(genpath('utils'));
addpath('constants');
addpath('initial_conditions');

nodes_file = FilePaths.NODES;
cells_file = FilePaths.CELLS;
bc_files = FilePaths.BOUNDARY_CONDITIONS;

cells = mesh_processor(nodes_file, cells_file, bc_files);

centroids = reshape([cells.centroid], 2, [])';
centroids_x = centroids(:, 1);

w0 = Config.INITIAL_CONDITIONS(centroids_x);

num_cells = length(cells);

rho0 = w0(1:num_cells);
rhou0 = w0(num_cells + 1:2 * num_cells);
E0 = w0(2 * num_cells + 1:3 * num_cells);

u0 = rhou0 ./ rho0;
p0 = (Air.GAMMA - 1) .* (E0 - 0.5 .* (rhou0 .^ 2) ./ rho0);

fprintf('num_cells      = %d\n', num_cells);
fprintf('rho0 min/max   = %.6f / %.6f\n', min(rho0), max(rho0));
fprintf('u0   min/max   = %.6f / %.6f\n', min(u0), max(u0));
fprintf('p0   min/max   = %.6f / %.6f\n', min(p0), max(p0));
fprintf('E0   min/max   = %.6f / %.6f\n', min(E0), max(E0));

figure;

subplot(3, 1, 1);
plot(centroids_x, rho0, 'o-');
grid on;
xlabel('x');
ylabel('\rho');
title('Initial density');

subplot(3, 1, 2);
plot(centroids_x, u0, 'o-');
grid on;
xlabel('x');
ylabel('u');
title('Initial velocity');

subplot(3, 1, 3);
plot(centroids_x, p0, 'o-');
grid on;
xlabel('x');
ylabel('p');
title('Initial pressure');
