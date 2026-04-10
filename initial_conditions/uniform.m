function state = uniform(pressure, density, velocity, cells_centroid_x)
    momentum = density * velocity;
    energy = get_internal_energy(density, pressure, velocity);

    num_cells = size(cells_centroid_x, 1);
    pressure_vec = ones(1, num_cells) * pressure;
    momentum_vec = ones(1, num_cells) * momentum;
    energy_vec = ones(1, num_cells) * energy;

    state = [pressure_vec; momentum_vec; energy_vec];
end
