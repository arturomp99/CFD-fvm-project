clc; clear all; close all;
%% Configuration

addpath('mesh_processing');
addpath('utils');
addpath('constants');
addpath('initial_conditions');

nodes_file = FilePaths.NODES;
cells_file = FilePaths.CELLS;
bc_files = FilePaths.BOUNDARY_CONDITIONS;

%% Mesh processing

[mesh] = mesh_processor(nodes_file, cells_file, bc_files);

%% Problem configuration
% se concatenarán 3 vectores uno detrás de otro
%   - la densidad media para cada celda
%   - la cantidad de movimiento para cada celda
%   - la energía interna más la energía cinética para cada celda.

state_vector = zeros(3, mesh.num_cells);
initial_conditions = Config.INITIAL_CONDITIONS;
source_terms = Config.SOURCE_TERMS;
boundary_conditions = Config.BOUNDARY_CONDITIONS;
