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

    plot_time_variations(centroids_x, results);

    plot_at_specific_time(centroids_x, results, 0.2);

end
