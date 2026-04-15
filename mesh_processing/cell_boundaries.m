function [boundary_faces, boundary_surface_ids] = cell_boundaries(cell_node_indices, bcs_data)
    %CELL_BOUNDARIES Identify which faces of a cell lie on boundary surfaces.
    %   [boundary_faces, boundary_surface_ids] = CELL_BOUNDARIES(cell_node_indices, bcs_data)
    %   determines which faces of a cell belong to boundary surfaces by checking
    %   if both nodes of each face are members of a boundary surface.
    %
    %   Inputs:
    %   -------
    %   cell_node_indices : vector
    %       Node indices defining the cell polygon (in order).
    %   bcs_data : cell array
    %       Cell array where bcs_data{i} contains the node indices belonging
    %       to boundary surface i.
    %
    %   Outputs:
    %   --------
    %   boundary_faces : vector
    %       Indices of faces that lie on boundary surfaces.
    %   boundary_surface_ids : vector
    %       For each boundary face, the index of the boundary surface it belongs to.
    %       Same length as boundary_faces.
    %
    %   Example:
    %   --------
    %   If cell has 4 faces and face 2 is on boundary surface 3:
    %       boundary_faces = [2]
    %       boundary_surface_ids = [3]

    num_faces = length(cell_node_indices);
    num_bc_surfaces = length(bcs_data);
    
    boundary_faces = [];
    boundary_surface_ids = [];
    
    for face_idx = 1:num_faces
        % Get the two node indices that form this face
        [node1_idx, node2_idx] = get_face_node_indices(cell_node_indices, face_idx);
        
        % Check if both nodes belong to any boundary surface
        for bc_idx = 1:num_bc_surfaces
            bc_nodes = bcs_data{bc_idx};
            if is_face_on_boundary(node1_idx, node2_idx, bc_nodes)
                boundary_faces = [boundary_faces, face_idx];
                boundary_surface_ids = [boundary_surface_ids, bc_idx];
                break;  % Face can only belong to one BC surface
            end
        end
    end
end
