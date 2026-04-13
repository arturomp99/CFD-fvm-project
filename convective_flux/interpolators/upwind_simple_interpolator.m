function [A, b] = upwind_simple_interpolator(state, cells)
    gamma = Air.GAMMA;
    num_cells = length(cells);
    densities = state(1:N);
    momentums = state(N + 1:2 * N);
    energies = state(2 * N + 1:3 * N);
    velocities = momentums ./ densities;
    pressures = (gamma - 1) * (energies - 0.5 .* densities .* velocities ^ 2);
    enthalpies = (energies + pressures) ./ densities;

    % initialize the matrices
    A = sparse(3 * N, 3 * N);
    b = zeros(3 * N, 1);

    for i = 1:num_cells
        % Find the Jacobian
        cell_density = densities(i);
        cell_momentum = momentums(i);
        cell_energy = energies(i);
        cell_velocity = velocities(i);
        cell_pressure = pressures(i);
        cell_enthalpy = enthalpies(i);

        cell_jacobian = [
                         0, 1, 0;
                         (gamma - 3) / 2 * cell_velocity ^ 2, (3 - gamma) * cell_velocity, gamma - 1;
                         (gamma - 1) / 2 * cell_velocity ^ 3 - cell_velocity * cell_enthalpy, cell_enthalpy - (gamma - 1) * cell_velocity ^ 2, gamma * cell_velocity
                         ];

        assembly_indices = (i - 1) * 3 + (1:3);

        % Bloque diagonal
        A(assembly_indices, assembly_indices) = ...
            A(assembly_indices, assembly_indices) - cell_jacobian / dx;

        % Bloque subdiagonal

    end

end
