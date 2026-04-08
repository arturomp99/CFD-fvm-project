function energy = get_internal_energy(density, pressure, velocity)
    energy = density * (Air.C_V * get_temperature(density, pressure) + velocity ^ 2/2);
end
