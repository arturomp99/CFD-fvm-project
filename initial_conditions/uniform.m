function state = uniform(pressure, density, velocity, cells_centroid_x)
    %UNIFORM Generates uniform initial conditions across all cells.
    %   state = UNIFORM(pressure, density, velocity, cells_centroid_x) creates a
    %   constant initial state vector [density; momentum; total energy] with the
    %   same primitive values in every cell.
    %
    %   Inputs:
    %   -------
    %   pressure         : double - Uniform static pressure. [Pa]
    %   density          : double - Uniform mass density. [kg/m^3]
    %   velocity         : double - Uniform flow velocity. [m/s]
    %   cells_centroid_x : column vector (N x 1) - x-coordinates of cell centroids (used to determine N). [m]
    %
    %   Outputs:
    %   --------
    %   state : column vector (3*N x 1)
    %     Concatenated initial state [density (N); momentum (N); total energy (N)].

    momentum = density * velocity;
    energy = get_internal_energy(density, pressure, velocity);

    num_cells = size(cells_centroid_x, 1);
    density_vec = ones(num_cells, 1) * density;
    momentum_vec = ones(num_cells, 1) * momentum;
    energy_vec = ones(num_cells, 1) * energy;

    state = [density_vec; momentum_vec; energy_vec];
end
