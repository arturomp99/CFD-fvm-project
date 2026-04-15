function U_bc = apply_bc_state(U_interior, bc_type, normal_dir, bc_params)
    %APPLY_BC_STATE Compute boundary/ghost state based on BC type.
    %   U_bc = APPLY_BC_STATE(U_interior, bc_type, normal_dir, bc_params) returns the
    %   ghost cell state for applying boundary conditions in flux calculations.
    %
    %   Inputs:
    %   -------
    %   U_interior : column vector [rho; rhou; E]
    %       State vector of the interior cell adjacent to the boundary.
    %   bc_type : string
    %       Boundary condition type:
    %         'open'     - Transmissive/open boundary (waves pass through)
    %         'wall'     - Solid wall boundary (reflective, zero normal velocity)
    %         'velocity' - Constant velocity boundary (fixed velocity imposed)
    %   normal_dir : scalar (+1 or -1)
    %       Direction of outward normal:
    %         +1 if boundary is to the right of the cell
    %         -1 if boundary is to the left of the cell
    %   bc_params : struct (optional)
    %       Additional parameters for specific BC types:
    %         .velocity - velocity value for 'velocity' BC [m/s]
    %
    %   Outputs:
    %   --------
    %   U_bc : column vector [rho; rhou; E]
    %       Boundary/ghost state for flux computation.
    %
    %   Notes:
    %   ------
    %   - 'open': Extrapolates the interior state (non-reflecting BC).
    %   - 'wall': Reflects momentum to enforce zero normal velocity.
    %   - 'velocity': Sets momentum to rho * u_bc for prescribed velocity.
    
    if nargin < 4
        bc_params = struct();
    end
    
    rho = U_interior(1);
    rhou = U_interior(2);
    E = U_interior(3);
    
    gamma = Air.GAMMA;
    
    switch bc_type
        case 'open'
            % Transmissive: extrapolate interior state (waves pass through)
            U_bc = U_interior;
            
        case 'wall'
            % Solid wall: reflect momentum (zero normal velocity)
            % Mirror the velocity component normal to the wall
            U_bc = [rho; -rhou; E];
            
        case 'velocity'
            % Constant velocity: impose fixed velocity at boundary
            % Extrapolate density and pressure from interior,
            % but set velocity to the prescribed value
            if ~isfield(bc_params, 'velocity')
                error('Velocity BC requires bc_params.velocity to be specified');
            end
            
            u_bc = bc_params.velocity;
            
            % Extrapolate density from interior
            rho_bc = rho;
            
            % Compute momentum with prescribed velocity
            rhou_bc = rho_bc * u_bc;
            
            % Extrapolate pressure from interior, then compute energy
            u_int = rhou / max(rho, 1e-10);
            p_int = (gamma - 1) * (E - 0.5 * rho * u_int^2);
            p_int = max(p_int, 1e-10);
            
            % Energy with prescribed velocity and extrapolated pressure
            E_bc = p_int / (gamma - 1) + 0.5 * rho_bc * u_bc^2;
            
            U_bc = [rho_bc; rhou_bc; E_bc];
            
        otherwise
            % Default to transmissive
            U_bc = U_interior;
    end
end
