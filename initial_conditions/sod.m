function state = sod(pressure, density, velocity, shock_pos, cells_centroid_x)
    %SOD Generates initial conditions for the Sod shock tube problem.
    %   state = SOD(pressure, density, velocity, shock_pos, cells_centroid_x)
    %   assigns left-state and right-state primitive variables on either side of
    %   a discontinuity at x = shock_pos, then converts them to the conserved
    %   variable vector [density; momentum; total energy].
    %
    %   Inputs:
    %   -------
    %   pressure : struct with fields .left and .right - Pressure values [Pa].
    %   density  : struct with fields .left and .right - Density values [kg/m^3].
    %   velocity : struct with fields .left and .right - Velocity values [m/s].
    %   shock_pos : double - x-coordinate of the initial discontinuity. [m]
    %   cells_centroid_x : column vector (N x 1) - x-coordinates of cell centroids. [m]
    %
    %   Outputs:
    %   --------
    %   state : column vector (3*N x 1)
    %     Concatenated initial state [density (N); momentum (N); total energy (N)].

    % Example usage in the Config.m file
    % INITIAL_CONDITIONS = @(pos) ... % funcion de la posición
    %         sod( ...
    %         struct('left', 1., 'right', 0.1), ... % pressure
    %         struct('left', 1., 'right', 0.125), ... % density
    %         struct('left', 0., 'right', 0.0), ... % velocity
    %         0.5, ...
    %         pos ...
    %     );

    momentum.left = density.left * velocity.left;
    momentum.right = density.right * velocity.right;

    energy.left = get_internal_energy(density.left, pressure.left, velocity.left);
    energy.right = get_internal_energy(density.right, pressure.right, velocity.right);

    num_cells = size(cells_centroid_x, 1);

    density_vec = zeros(num_cells, 1);
    momentum_vec = zeros(num_cells, 1);
    energy_vec = zeros(num_cells, 1);

    left_mask = (cells_centroid_x <= shock_pos);
    right_mask = (cells_centroid_x > shock_pos);

    density_vec(left_mask) = density.left;
    density_vec(right_mask) = density.right;

    momentum_vec(left_mask) = momentum.left;
    momentum_vec(right_mask) = momentum.right;

    energy_vec(left_mask) = energy.left;
    energy_vec(right_mask) = energy.right;

    state = [density_vec; momentum_vec; energy_vec];
end
