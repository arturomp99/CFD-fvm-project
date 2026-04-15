function [cells, boundary_info] = mesh_processor(nodes_file, cells_file, bc_files)
    %MESH_PROCESSOR Load and process mesh data from files.
    %
    %   Outputs:
    %   --------
    %   cells : struct array
    %       Cell structures with geometry, connectivity, and boundary info.
    %   boundary_info : struct
    %       Structure containing boundary surface information and BC types.
    
    % Loading the data...
    [nodes_data, cells_data, bcs_data] = mesh_data_loader( ...
        nodes_file, cells_file, bc_files);
    % Data loaded. Now, process the geometry.
    [cells, boundary_info] = process_mesh_geometry(nodes_data, cells_data, bcs_data);
    % Finally, compute the connectivity
    cells = compute_connectivity(cells);
    
    % Assign boundary condition configuration from Config
    boundary_info.boundary_types = Config.BOUNDARY_TYPES;
    boundary_info.boundary_velocities = Config.BOUNDARY_VELOCITIES;
end
