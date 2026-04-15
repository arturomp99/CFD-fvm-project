function [cells, boundary_info] = process_mesh_geometry( ...
        nodes_data, ...
        cells_data, ...
        bcs_data ...
    )
    %PROCESS_MESH_GEOMETRY Process mesh geometry and identify boundary faces.
    %
    %   Outputs:
    %   --------
    %   cells : struct array
    %       Cell structures with geometry and boundary face information.
    %   boundary_info : struct
    %       Structure containing:
    %         .surface_names - cell array of boundary surface identifiers
    %         .surface_nodes - cell array of node indices for each surface

    % initialize the cells vector
    num_cells = length(cells_data);
    cell_data = struct( ...
        'nodes', [], ...
        'faces', [], ...
        'volume', 0, ...
        'centroid', [0, 0], ...
        'area', [], ...
        'normals', [0, 0], ...
        'boundary_faces', [], ...  % indices of faces that are boundaries
        'boundary_surface_ids', [] ... % which boundary surface each boundary face belongs to
    );
    cells = repmat(cell_data, 1, num_cells);

    % Process boundary surface node sets
    num_bc_surfaces = length(bcs_data);
    boundary_info = struct();
    boundary_info.surface_names = cell(num_bc_surfaces, 1);
    boundary_info.surface_nodes = cell(num_bc_surfaces, 1);
    
    for bc_idx = 1:num_bc_surfaces
        boundary_info.surface_names{bc_idx} = sprintf('surface_%d', bc_idx);
        boundary_info.surface_nodes{bc_idx} = bcs_data{bc_idx};
    end

    for cell_index = 1:num_cells
        % Process the data of one single cell (including boundary faces)
        cells(cell_index) = process_cell(cells_data{cell_index}, nodes_data, bcs_data);
    end

end
