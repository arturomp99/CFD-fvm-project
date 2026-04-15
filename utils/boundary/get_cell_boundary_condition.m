function [bc_type, bc_params] = get_cell_boundary_condition(cell, side, boundary_info)
    %GET_CELL_BOUNDARY_CONDITION Find the BC type for a cell's boundary face.
    %   [bc_type, bc_params] = GET_CELL_BOUNDARY_CONDITION(cell, side, boundary_info)
    %   determines which boundary surface a cell's face belongs to and returns
    %   the corresponding boundary condition type and parameters.
    %
    %   Inputs:
    %   -------
    %   cell : struct
    %       Cell structure with .nodes (Vx2 coordinates) and .boundary_surface_ids
    %   side : string
    %       'left' or 'right' - which face to check (1D: lower or higher x)
    %   boundary_info : struct
    %       Structure with:
    %         .surface_nodes     - cell array of node indices for each surface
    %         .boundary_types    - cell array of BC types for each surface
    %         .boundary_velocities - array of velocity values for each surface
    %
    %   Outputs:
    %   --------
    %   bc_type : string
    %       'open', 'wall', or 'velocity'. Defaults to 'open'.
    %   bc_params : struct
    %       Additional parameters (e.g., .velocity for 'velocity' BC).
    
    % Default values
    bc_type = 'open';
    bc_params = struct();
    
    % Validate inputs
    if isempty(boundary_info) || ~isfield(boundary_info, 'boundary_types')
        return;
    end
    
    if isempty(boundary_info.boundary_types)
        return;
    end
    
    if isempty(cell.boundary_surface_ids)
        return;
    end
    
    % Find which boundary surface corresponds to the requested side
    bc_surface_id = find_boundary_surface_for_side(cell, side);
    
    % Get the BC type and parameters for that surface
    [bc_type, bc_params] = get_bc_from_surface_id(bc_surface_id, boundary_info);
end
