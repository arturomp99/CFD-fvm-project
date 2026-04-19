function cells = mesh_processor(nodes_file, cells_file, bc_files)
    %MESH_PROCESSOR Load and process mesh data from files.
    %
    %   Outputs:
    %   --------
    %   cells : struct array
    %       Cell structures with geometry, connectivity, and boundary info.

    % Loading the data...
    [nodes_data, cells_data, bcs_data] = mesh_data_loader( ...
        nodes_file, cells_file, bc_files);
    % Data loaded. Now, process the geometry.
    cells = process_mesh_geometry(nodes_data, cells_data, bcs_data);
    % Finally, compute the connectivity
    cells = compute_connectivity(cells);
end
