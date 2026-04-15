function [bc_type, bc_params] = get_bc_from_surface_id(bc_surface_id, boundary_info)
    %GET_BC_FROM_SURFACE_ID Get boundary condition type and parameters for a surface.
    %   [bc_type, bc_params] = GET_BC_FROM_SURFACE_ID(bc_surface_id, boundary_info)
    %   retrieves the boundary condition type and any additional parameters
    %   (like velocity) for the specified boundary surface.
    %
    %   Inputs:
    %   -------
    %   bc_surface_id : integer
    %       Index of the boundary surface in boundary_info.
    %   boundary_info : struct
    %       Structure with:
    %         .boundary_types    - cell array of BC types for each surface
    %         .boundary_velocities - array of velocity values for each surface
    %
    %   Outputs:
    %   --------
    %   bc_type : string
    %       'open', 'wall', or 'velocity'.
    %   bc_params : struct
    %       Additional parameters (e.g., .velocity for 'velocity' BC).

    bc_type = 'open';  % Default
    bc_params = struct();
    
    if isempty(bc_surface_id)
        return;
    end
    
    boundary_types = boundary_info.boundary_types;
    
    if bc_surface_id > length(boundary_types)
        return;
    end
    
    bc_type = boundary_types{bc_surface_id};
    
    % Get velocity parameter if this is a velocity BC
    if strcmp(bc_type, 'velocity')
        if isfield(boundary_info, 'boundary_velocities') && ...
           bc_surface_id <= length(boundary_info.boundary_velocities)
            bc_params.velocity = boundary_info.boundary_velocities(bc_surface_id);
        else
            error('Velocity BC on surface %d requires boundary_velocities', bc_surface_id);
        end
    end
end
