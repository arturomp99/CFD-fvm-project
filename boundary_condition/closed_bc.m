function boundary_state = closed_bc(face_state, gamma)
    % CLOSED_BC - Closed boundary condition (no flow through boundary)
    %
    %   Sets velocity = 0 at the boundary while preserving the thermodynamic
    %   state (density and pressure) from the interior. This represents a
    %   solid wall or closed boundary with no mass flow.
    %
    % Input
    %   face_state  - State vector from interior [ρ; ρu; E] [3×1]
    %   gamma       - Specific heat ratio γ = cp/cv [-]
    %
    % Output
    %   boundary_state - State with zero velocity [ρ; 0; E'] [3×1]

    % Extract interior state
    rho = face_state(1);
    rhou = face_state(2);
    E = face_state(3);

    % Current velocity
    u_interior = rhou / rho;

    % For closed boundary: velocity = 0, density = interior density
    u_boundary = 0.0;

    % Compute kinetic energy difference
    ke_interior = 0.5 * rho * u_interior ^ 2;
    ke_boundary = 0.5 * rho * u_boundary ^ 2;

    % Adjust total energy: E_boundary = E_interior - ke_interior + ke_boundary
    E_boundary = E - ke_interior + ke_boundary;

    % Return boundary state with zero velocity
    boundary_state = [rho; rho * u_boundary; E_boundary];
end
