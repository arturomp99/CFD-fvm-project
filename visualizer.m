function visualizer(results, centroids_x)
    %VISUALIZER Plots the space-time evolution of the 1D Euler solution.
    %   VISUALIZER(results, centroids_x) produces three space-time contour plots,
    %   one for each conserved variable: density, momentum, and total energy.
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

    density = results(:, 2:N + 1);
    momentum = results(:, N + 2:2 * N + 1);
    energy = results(:, 2 * N + 2:3 * N + 1);

    [X, T] = meshgrid(centroids_x, t_vec);

    figure;

    subplot(3, 1, 1);

    if (Config.DIFFUMINATED_VISUALIZATION)
        contourf(X, T, density, 20);
    else
        imagesc(centroids_x, t_vec, density);
    end

    colorbar();
    xlabel('x [m]');
    ylabel('t [s]');
    title('Density \rho [kg/m^3]');

    subplot(3, 1, 2);

    if (Config.DIFFUMINATED_VISUALIZATION)
        contourf(X, T, momentum, 20);
    else
        imagesc(centroids_x, t_vec, momentum);
    end

    colorbar();
    xlabel('x [m]');
    ylabel('t [s]');
    title('Momentum \rho v [kg/(m^2\cdots)]');

    subplot(3, 1, 3);

    if (Config.DIFFUMINATED_VISUALIZATION)
        contourf(X, T, energy, 20);
    else
        imagesc(centroids_x, t_vec, energy);
    end

    colorbar();
    xlabel('x [m]');
    ylabel('t [s]');
    title('Total energy E [J/m^3]');

end
