function state = sod(pressure, density, velocity, shock_pos, cells_centroid_x)
    % SOD shock tube problem initial conditions

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

    density_vec(cells_centroid_x < shock_pos) = density.left;
    density_vec(cells_centroid_x > shock_pos) = density.right;
    momentum_vec(cells_centroid_x < shock_pos) = momentum.left;
    momentum_vec(cells_centroid_x > shock_pos) = momentum.right;
    energy_vec(cells_centroid_x < shock_pos) = energy.left;
    energy_vec(cells_centroid_x > shock_pos) = energy.right;

    state = [density_vec; momentum_vec; energy_vec];
end
