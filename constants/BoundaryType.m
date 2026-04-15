classdef BoundaryType
    %BOUNDARYTYPE Enumeration of supported boundary condition types.
    %
    %   Available types:
    %   ----------------
    %   OPEN     - Transmissive/open boundary. Waves pass through without
    %              reflection. Uses extrapolation from interior.
    %   WALL     - Solid wall boundary. Reflective condition with zero
    %              normal velocity. Momentum is reflected, density and
    %              energy are extrapolated.
    %   VELOCITY - Constant velocity boundary. Imposes a fixed velocity
    %              at the boundary. Requires specifying the velocity value
    %              in Config.BOUNDARY_VELOCITIES.
    
    properties (Constant)
        OPEN = 'open';
        WALL = 'wall';
        VELOCITY = 'velocity';
    end
    
    methods (Static)
        function valid = isValid(type_str)
            %ISVALID Check if a string is a valid boundary type.
            valid = strcmp(type_str, BoundaryType.OPEN) || ...
                    strcmp(type_str, BoundaryType.WALL) || ...
                    strcmp(type_str, BoundaryType.VELOCITY);
        end
        
        function U_bc = applyBoundaryState(U_interior, bc_type, normal_dir)
            %APPLYBOUNDARYSTATE Compute ghost/boundary state based on BC type.
            %
            %   Inputs:
            %   -------
            %   U_interior : [rho; rhou; E] state vector of interior cell
            %   bc_type    : string, either 'open' or 'wall'
            %   normal_dir : +1 or -1, direction of outward normal
            %                (+1 if boundary is to the right, -1 if to the left)
            %
            %   Outputs:
            %   --------
            %   U_bc : [rho; rhou; E] boundary/ghost state
            
            rho = U_interior(1);
            rhou = U_interior(2);
            E = U_interior(3);
            
            switch bc_type
                case BoundaryType.OPEN
                    % Transmissive: extrapolate interior state
                    U_bc = U_interior;
                    
                case BoundaryType.WALL
                    % Solid wall: reflect momentum, keep density and energy
                    % For 1D: reverse the momentum to simulate reflection
                    U_bc = [rho; -rhou; E];
                    
                case BoundaryType.VELOCITY
                    % Velocity BC requires additional parameters
                    % Use apply_bc_state utility function instead
                    error('Use apply_bc_state with bc_params for velocity BC');
                    
                otherwise
                    error('Unknown boundary condition type: %s', bc_type);
            end
        end
    end
end
