clc; clear all; close all;
%% Configuration

addpath('mesh_processing');
addpath(genpath('utils'));
addpath('constants');
addpath('initial_conditions');
addpath('problems/fvm_1D_euler');
addpath(genpath('convective_flux'));
addpath('stopping_criteria')
addpath('timestep_control')
addpath('propagators')
 addpath('results_manager')

nodes_file = FilePaths.NODES;
cells_file = FilePaths.CELLS;
bc_files = FilePaths.BOUNDARY_CONDITIONS;

%% Mesh processing

cells = mesh_processor(nodes_file, cells_file, bc_files);

%% Problem configuration
% se concatenarán 3 vectores uno detrás de otro
%   - la densidad media para cada celda
%   - la cantidad de movimiento para cada celda
%   - la energía interna más la energía cinética para cada celda.

initial_conditions = Config.INITIAL_CONDITIONS;
centroids_x = reshape([cells.centroid], 2, [])';
centroids_x = centroids_x(:, 1);
w0 = initial_conditions(centroids_x);

problem = @(state, time) fvm_1D_euler(state, cells);

time_integrator = @bw_euler;

timestep_calculator = @(w, t) constant_dt(w, t, 0.01);

stopping_condition = @(w, t) stop_at_time(t, 1.5);

% Use a sampling period of 0.01 s. Configure the results manager using an
% anonymous function, so the solver only sees a function that depends on
% the state vector, the time and the previous results structure.
manager = @(w, t, old_results) sample_results(w, t, old_results, 0.01);

%% Resolution
results = solver( ...
    w0, ...
    0., ...
    problem, ...
    time_integrator, ...
    timestep_calculator, ...
    stopping_condition, ...
    manager ...
);
