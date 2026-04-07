function [mesh_data] = mesh_processor(nodes_file, cells_file, bc_files)
    % Loading the data...
    [nodes_data, cells_data, bcs_data] = mesh_data_loader( ...
        nodes_file, cells_file, bc_files);
    % Data loaded. Now, process the geometry.
    mesh_data = process_mesh_geometry(nodes_data, cells_data, bcs_data);
    % Finally, compute the connectivity
    % TODO: Implement this ----
    % mesh_data.connectivity = compute_connectivity(mesh_data);
    % TODO: ----
end
