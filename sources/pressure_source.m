function [A, b] = pressure_source(state, centroids_x, t, source_position, pressure_forcing, source_width)
    %PRESSURE_SOURCE Applies a localised pressure forcing at a given x position.
    %
    %   [A, b] = PRESSURE_SOURCE(state, centroids_x, t, source_position, pressure_forcing, source_width)
    %   returns source terms that add pressure forcing to cells near the specified position.
    %   This models a pressure disturbance or actuator that pushes fluid at a fixed point in space.
    %
    %   The forcing affects the momentum equation directly and adjusts energy for consistency.
    %
    %   Inputs:
    %   -------
    %   state           : column vector (3*N x 1) - Current conserved state
    %   centroids_x     : column vector (N x 1)   - x-coordinates of cell centroids [m]
    %   t               : double                  - Current time (can be used for time-varying forcing) [s]
    %   source_position : double                  - x-coordinate where pressure forcing is applied [m]
    %   pressure_forcing: double                  - Pressure forcing magnitude [Pa]
    %   source_width    : double                  - Width of the forcing region [m]
    %
    %   Outputs:
    %   --------
    %   A : sparse matrix (3N x 3N) - Matrix contribution (zero for explicit sources)
    %   b : column vector (3N x 1)  - Source term vector

    num_cells = length(centroids_x);
    source_vec = zeros(num_cells * 3, 1);
    [A, b] = initialize_A_b(num_cells);

    % Find cells within source_width of the source position
    distance_from_source = abs(centroids_x - source_position);
    source_cells = distance_from_source <= (source_width / 2);

    if ~any(source_cells)
        return; % No cells in source region
    end

    % Extract current state for source cells
    [rho, rhou, E] = state_vec2states(state);
    rho_source = rho(source_cells);

    % Pressure forcing creates a momentum source: dp/dx ≈ pressure_forcing / source_width
    % This is a simplified model - in reality pressure gradients drive flow
    momentum_forcing = pressure_forcing * ones(sum(source_cells), 1);

    % Apply forcing to momentum equation (pressure gradient drives momentum)
    momentum_indices = num_cells + find(source_cells);
    source_vec(momentum_indices) = momentum_forcing;

    % For thermodynamic consistency, adjust energy equation
    % Pressure work: dE/dt includes -p * div(u) + work terms
    % Simplified: add pressure forcing * velocity to energy equation
    u_source = rhou(source_cells) ./ rho_source;
    energy_forcing = pressure_forcing * u_source;

    energy_indices = 2 * num_cells + find(source_cells);
    source_vec(energy_indices) = energy_forcing;

    if Config.IS_SOURCE_IMPLICIT
        A = diag(source_vec);
    else
        b = source_vec;
    end

end
