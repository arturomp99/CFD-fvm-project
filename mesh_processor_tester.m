%% Configuration

addpath('mesh_processing');
addpath('utils');

nodes_file = 'L_grid/nodes_L_1.dat';
cells_file = 'L_grid/cells_L_1.dat';
bc_files = {'L_grid/bc_L_1_1.dat', 'L_grid/bc_L_1_2.dat', ...
                'L_grid/bc_L_1_3.dat', 'L_grid/bc_L_1_4.dat', ...
                'L_grid/bc_L_1_5.dat', 'L_grid/bc_L_1_6.dat', ...
            'L_grid/bc_L_1_7.dat'};

%% Mesh processing

[mesh] = mesh_processor(nodes_file, cells_file, bc_files);
