function [A, b] = point_velocity_source(centroids_x, source_position, source_strength)
    %POINT_VELOCITY_SOURCE Applies a localised momentum forcing at a given x position.
    %   source = POINT_VELOCITY_SOURCE(state, centroids_x, t, source_position, source_strength)
    %   returns a source term vector that adds momentum (a body force) to the
    %   single cell closest to source_position. The density and energy equations
    %   are unforced (mass and heat are not injected, only a force is applied).
    %
    %   This models a device such as a fan or actuator disk that pushes the fluid
    %   at a fixed point in space.
    %
    %   Inputs:
    %   -------
    %   state           : column vector (3*N x 1) - Current conserved state (unused here).
    %   centroids_x     : column vector (N x 1)   - x-coordinates of cell centroids. [m]
    %   t               : double                  - Current time (can be used for time-varying forcing). [s]
    %   source_position : double                  - x-coordinate where forcing is applied. [m]
    %   source_strength : double                  - Momentum forcing magnitude. [kg/(m^2 s^2)]
    %
    %   Outputs:
    %   --------
    %   source : column vector (3*N x 1)
    %     Source term vector. Only the momentum component of the nearest cell
    %     is non-zero.

    num_cells = length(centroids_x);
    [A, b] = initialize_A_b(num_cells);

    % Find the cell closest to the requested source position
    [~, source_cell] = min(abs(centroids_x - source_position));

    % Add the forcing to the momentum equation of that cell
    b(num_cells + source_cell) = source_strength;
end
