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
        % Process the data of one single cell and save it in mesh_data.
        cells(cell_index) = process_cell(cells_data{cell_index}, nodes_data);
        
        % Identify boundary faces for this cell
        cell_node_indices = cells_data{cell_index};
        num_faces = length(cell_node_indices);
        
        boundary_faces = [];
        boundary_surface_ids = [];
        
        for face_idx = 1:num_faces
            % Get the two node indices that form this face
            node1_idx = cell_node_indices(face_idx);
            node2_idx = cell_node_indices(mod(face_idx, num_faces) + 1);
            
            % Check if both nodes belong to any boundary surface
            for bc_idx = 1:num_bc_surfaces
                bc_nodes = bcs_data{bc_idx};
                if ismember(node1_idx, bc_nodes) && ismember(node2_idx, bc_nodes)
                    boundary_faces = [boundary_faces, face_idx]; %#ok<AGROW>
                    boundary_surface_ids = [boundary_surface_ids, bc_idx]; %#ok<AGROW>
                    break;  % Face can only belong to one BC surface
                end
            end
        end
        
        cells(cell_index).boundary_faces = boundary_faces;
        cells(cell_index).boundary_surface_ids = boundary_surface_ids;
    end

end
