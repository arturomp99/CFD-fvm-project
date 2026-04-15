function bc_surface_id = find_boundary_surface_for_side(cell, side)
    %FIND_BOUNDARY_SURFACE_FOR_SIDE Find which boundary surface corresponds to a side.
    %   bc_surface_id = FIND_BOUNDARY_SURFACE_FOR_SIDE(cell, side) determines
    %   which boundary surface the left or right face of the cell belongs to
    %   by checking x-coordinates.
    %
    %   Inputs:
    %   -------
    %   cell : struct
    %       Cell structure with .nodes, .boundary_faces, and .boundary_surface_ids
    %   side : string
    %       'left' or 'right' - which side to check (lower or higher x)
    %
    %   Outputs:
    %   --------
    %   bc_surface_id : integer or empty
    %       Index of the boundary surface, or empty if no boundary face on that side.

    bc_surface_id = [];
    
    if isempty(cell.boundary_faces)
        return;
    end
    
    cell_nodes = cell.nodes;
    
    % Initialize target x based on side
    if strcmp(side, 'left')
        target_x = inf;   % Looking for minimum x
    else
        target_x = -inf;  % Looking for maximum x
    end
    
    % Search through boundary faces to find the one on the specified side
    for k = 1:length(cell.boundary_faces)
        face_idx = cell.boundary_faces(k);
        face_avg_x = get_face_average_x(cell_nodes, face_idx);
        
        is_better_match = (strcmp(side, 'left') && face_avg_x < target_x) || ...
                          (strcmp(side, 'right') && face_avg_x > target_x);
        
        if is_better_match
            target_x = face_avg_x;
            bc_surface_id = cell.boundary_surface_ids(k);
        end
    end
end
