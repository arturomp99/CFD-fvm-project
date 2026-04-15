function [bc_type, bc_params] = get_boundary_type(cell, side, boundary_info)
    %GET_BOUNDARY_TYPE Determine the boundary type for a cell's boundary face.
    %   [bc_type, bc_params] = GET_BOUNDARY_TYPE(cell, side, boundary_info) returns the
    %   boundary condition type ('open', 'wall', or 'velocity') and any additional
    %   parameters for the specified cell.
    %
    %   Inputs:
    %   -------
    %   cell : struct
    %       Cell structure with .boundary_surface_ids field.
    %   side : string
    %       'left' or 'right' indicating which face to check.
    %   boundary_info : struct
    %       Structure with:
    %         .boundary_types - cell array specifying BC type for each surface
    %         .boundary_velocities - (optional) array of velocity values for 
    %                                'velocity' BC [m/s]
    %
    %   Outputs:
    %   --------
    %   bc_type : string
    %       'open' (transmissive), 'wall' (reflective), or 'velocity'.
    %       Defaults to 'open'.
    %   bc_params : struct
    %       Additional parameters for the BC (e.g., .velocity for 'velocity' BC).
    
    bc_type = 'open';  % Default: transmissive boundary
    bc_params = struct();
    
    if isempty(boundary_info) || ~isfield(boundary_info, 'boundary_types')
        return;
    end
    
    boundary_types = boundary_info.boundary_types;
    if isempty(boundary_types)
        return;
    end
    
    % Check if this cell has boundary faces and get the appropriate BC type
    if ~isempty(cell.boundary_surface_ids)
        % For 1D, we check if this is a left or right boundary
        % Use the first boundary surface ID found for this cell
        bc_surface_id = cell.boundary_surface_ids(1);
        if bc_surface_id <= length(boundary_types)
            bc_type = boundary_types{bc_surface_id};
            
            % Get velocity parameter if this is a velocity BC
            if strcmp(bc_type, 'velocity')
                if isfield(boundary_info, 'boundary_velocities') && ...
                   bc_surface_id <= length(boundary_info.boundary_velocities)
                    bc_params.velocity = boundary_info.boundary_velocities(bc_surface_id);
                else
                    error('Velocity BC on surface %d requires boundary_velocities to be defined', bc_surface_id);
                end
            end
        end
    end
end
