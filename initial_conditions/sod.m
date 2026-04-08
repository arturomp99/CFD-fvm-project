function state = sod(pressure, density, velocity, shock_pos, cells_centroid_x)
    % SOD shock tube problem initial conditions

    momentum.left = pressure.left * velocity.left;
    momentum.right = pressure.right * velocity.right;
    energy.left = get_internal_energy(density.left, pressure.left, velocity.left);
    energy.right = get_internal_energy(density.right, pressure.right, velocity.right);

    num_cells = size(cells_centroid_x, 1);
    density_vec = zeros(num_cells, 1);
    momentum_vec = zeros(num_cells, 1);
    energy_vec = zeros(num_cells, 1);

    density_vec(cells_centroid_x < shock_pos) = density.left;
    density_vec(cells_centroid_x > shock_pos) = density.right;
    momentum_vec(cells_centroid_x < shock_pos) = momentum.left;
    momentum_vec(cells_centroid_x < shock_pos) = momentum.left;
    energy_vec(cells_centroid_x < shock_pos) = energy.left;
    energy_vec(cells_centroid_x < shock_pos) = energy.left;

    state = [density_vec; momentum_vec; energy_vec];
end
