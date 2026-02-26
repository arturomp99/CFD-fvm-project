function [nodes_data, cells_data, bcs_data] = mesh_data_loader(...
    nodes_file, cells_file, bc_files)
% MESH_DATA_LOADER Loads mesh data from mesh data files.
%
%   Inputs
%   ------
%   nodes_file : string
%     Name of the file that contains the coordinates of the nodes.
%   cells_data : string
%     Name of the file that contains the nodes forming each cell.
%   bcs_data : cell(string)
%     List of names of files that contain the nodes forming each boundary
%     condition surface.
%
%   Outputs
%   -------
%   nodes_data : cell(double)
%     Cell array with the coordinates of each node. The coordinates are
%     expressed in m.
%   cells_data : cell(double)
%     Cell array with the node indices forming each cell.
%   bcs_data : cell(double)
%     Cell array with the nodes forming each boundary condition surface.

    nodes_data = data_loader(nodes_file, "\t", 1, 1);
    cells_data = data_loader(cells_file, "\t", 1, 1);
    bcs_data = {1, length(bc_files)};
    for i=1:length(bc_files)
        current_bc_data = data_loader(bc_files{i}, "\t", 0, 1);
        bcs_data{i} = current_bc_data{1};
    end
end

