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

% cells is a (1 x N) struct array, where N is the number of mesh cells.
% Each element represents one cell and has the following fields:
%   .nodes        - (V x 2) array of node coordinates [x, y] for the V vertices of the cell. [m]
%   .faces        - (V x 2 x 2) array of face definitions; faces(i,1,:) and faces(i,2,:)
%                   are the [x,y] coordinates of the two endpoints of face i.
%   .volume       - scalar area of the cell (polyarea in 2D). [m^2]
%   .centroid     - (1 x 2) vector [cx, cy] of the geometric centroid. [m]
%   .area         - (1 x V) vector of edge lengths (face areas in 2D). [m]
%   .normals      - (V x 2) array of outward unit normal vectors for each face.
%   .connectivity - row vector of indices of the neighbouring cells that share a face.

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

propagator = Config.PROPAGATOR;

timestep_calculator = Config.TIMESTEP_CALCULATOR;

stopping_condition = Config.STOPPING_CONDITION;

% Use a sampling period of 0.01 s. Configure the results manager using an
% anonymous function, so the solver only sees a function that depends on
% the state vector, the time and the previous results structure.
manager = Config.RESULTS_MANAGER;

%% Resolution
results = solver( ...
    w0, ...
    0., ...
    problem, ...
    propagator, ...
    timestep_calculator, ...
    stopping_condition, ...
    manager ...
);

% results is a matrix of size (num_samples x (3*N + 1)), where num_samples
% is the number of recorded time steps and N is the number of mesh cells.
%
% Each row corresponds to one sampled instant:
%   results(i, 1)         -> time t [s]
%   results(i, 2:N+1)     -> density     for each cell [kg/m^3]
%   results(i, N+2:2N+1)  -> momentum    for each cell [kg/(m^2·s)]
%   results(i, 2N+2:3N+1) -> total energy for each cell [J/m^3]

%% Visualisation
visualizer(results, centroids_x);
