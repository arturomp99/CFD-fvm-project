clc; clear; close all;

addpath('mesh_processing');
addpath(genpath('utils'));
addpath('constants');

nodes_file = FilePaths.NODES;
cells_file = FilePaths.CELLS;
bc_files = FilePaths.BOUNDARY_CONDITIONS;

cells = mesh_processor(nodes_file, cells_file, bc_files);

num_cells = length(cells);
centroids = reshape([cells.centroid], 2, [])';

fprintf('Loaded cells: %d\n', num_cells);
fprintf('Centroid x min/max: %.6f / %.6f\n', min(centroids(:, 1)), max(centroids(:, 1)));
fprintf('Centroid y min/max: %.6f / %.6f\n', min(centroids(:, 2)), max(centroids(:, 2)));

figure;
plot(centroids(:, 1), centroids(:, 2), 'o');
grid on;
xlabel('x');
ylabel('y');
title('Cell centroids');
