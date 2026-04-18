function [A, b] = velocity_source_bc( ...
        state, ...
        centroids_x, ...
        t, ...
        target_velocity, ...
        relaxation_factor, ...
        boundary_width ...
    )
    %VELOCITY_SOURCE_BC Velocity source term near boundary
    %
    %   Applies a relaxation forcing to impose target velocity near the boundary.
    %   This acts as a "sponge layer" that smoothly forces the velocity towards
    %   the target value in cells near the boundary.
    %
    %   Input
    %   -----
    %   state            - Current state vector [3*N×1]
    %   centroids_x      - Cell centroid x-coordinates [N×1]
    %   t                - Current time [s] (unused)
    %   target_velocity  - Target velocity to impose [m/s]
    %   relaxation_factor- Relaxation strength (0-1, higher = stronger forcing)
    %   boundary_width   - Width of the forcing zone from left boundary [m]
    %
    %   Output
    %   ------
    %   A                - Matrix contribution (sparse, 3N×3N, typically zero)
    %   b                - Source term vector [3*N×1]

    num_cells = length(centroids_x);
    [A, b] = initialize_A_b(num_cells);
    source_vec = zeros(num_cells * 3, 1);

    % Find cells within boundary_width from the left boundary (x_min)
    x_min = min(centroids_x);
    boundary_cells = centroids_x <= (x_min + boundary_width);

    if ~any(boundary_cells)
        return; % No cells in boundary region
    end

    % Extract current state for all cells
    [rho, rhou, ~] = state_vec2states(state);

    % Current velocities in boundary cells
    u_current = rhou(boundary_cells) ./ rho(boundary_cells);

    % Target momentum for boundary cells
    rho_boundary = rho(boundary_cells);
    rhou_target = rho_boundary * target_velocity;

    % Relaxation forcing: source = relaxation_factor * (target - current)
    momentum_forcing = relaxation_factor * (rhou_target - rhou(boundary_cells));

    % Apply forcing to momentum equation
    source_indices = num_cells + find(boundary_cells);
    source_vec(source_indices) = momentum_forcing;

    % Also adjust energy to maintain consistency
    % The energy adjustment ensures thermodynamic consistency
    u_current_boundary = u_current;
    u_target = target_velocity;

    % Energy source = relaxation_factor * rho * 0.5 * (u_target^2 - u_current^2)
    energy_forcing = relaxation_factor * rho_boundary .* 0.5 .* (u_target .^ 2 - u_current_boundary .^ 2);

    energy_indices = 2 * num_cells + find(boundary_cells);
    source_vec(energy_indices) = energy_forcing;

    if Config.IS_SOURCE_IMPLICIT
        A = diag(source_vec);
    else
        b = source_vec;
    end

end
