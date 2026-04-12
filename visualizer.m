function visualizer(results, centroids_x)
    %VISUALIZER Plots the space-time evolution of the 1D Euler solution.
    %   VISUALIZER(results, centroids_x) produces a 3x2 grid of space-time plots,
    %   one for each of: density, velocity, momentum, temperature, pressure and
    %   total energy.
    %
    %   The horizontal axis is the cell centroid x-coordinate, the vertical axis
    %   is time, and the colour represents the variable's value at that (x, t) point.
    %
    %   Inputs:
    %   -------
    %   results : matrix (num_samples x (3*N + 1))
    %     Output of the solver as recorded by sample_results.
    %     Column 1 is time; the next three blocks of N columns are density,
    %     momentum, and total energy for all N cells.
    %   centroids_x : column vector (N x 1)
    %     x-coordinates of the cell centroids. [m]

    N = length(centroids_x);
    t_vec = results(:, 1);

    density  = results(:, 2 : N + 1);
    momentum = results(:, N + 2 : 2*N + 1);
    energy   = results(:, 2*N + 2 : 3*N + 1);

    velocity = momentum ./ density;
    internal_energy = energy - 0.5 * momentum .* velocity;  % e = E - 0.5*rho*v^2 [J/m^3]
    temperature = internal_energy ./ (density * Air.C_V);   % T = e/(rho*Cv) [K]
    pressure = density * Air.R .* temperature;              % p = rho*R*T [Pa]

    [X, T] = meshgrid(centroids_x, t_vec);

    vars   = {density, velocity, momentum, temperature, pressure, energy};
    titles = {'Density \rho [kg/m^3]', 'Velocity v [m/s]', ...
              'Momentum \rho v [kg/(m^2\cdots)]', 'Temperature T [K]', ...
              'Pressure p [Pa]', 'Total energy E [J/m^3]'};

    figure;

    for k = 1:6
        subplot(3, 2, k);

        if (Config.DIFFUMINATED_VISUALIZATION)
            contourf(X, T, vars{k}, 20);
        else
            imagesc(centroids_x, t_vec, vars{k});
            set(gca, 'YDir', 'normal');
        end

        colorbar();
        xlabel('x [m]');
        ylabel('t [s]');
        title(titles{k});
    end

end
