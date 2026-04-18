function plot_at_specific_time(centroids_x, results, t_target)
    centroids_x = centroids_x(:);
    N = length(centroids_x);

    t_vec = results(:, 1);

    density = results(:, 2:N + 1);
    momentum = results(:, N + 2:2 * N + 1);
    energy = results(:, 2 * N + 2:3 * N + 1);

    [~, idx_t] = min(abs(t_vec - t_target));
    t_plot = t_vec(idx_t);

    rho = density(idx_t, :).';
    rhou = momentum(idx_t, :).';
    E = energy(idx_t, :).';

    u = rhou ./ rho;
    p = (Air.GAMMA - 1) .* (E - 0.5 .* (rhou .^ 2) ./ rho);

    [x_sorted, idx] = sort(centroids_x);
    rho = rho(idx);
    u = u(idx);
    p = p(idx);

    figure;

    subplot(3, 1, 1);
    stairs(x_sorted, rho, 'b', 'LineWidth', 2);
    grid on;
    xlabel('x [m]');
    ylabel('\rho [kg/m^3]');
    title(sprintf('Density at t = %.4f s', t_plot));

    subplot(3, 1, 2);
    plot(x_sorted, u, 'b', 'LineWidth', 2);
    grid on;
    xlabel('x [m]');
    ylabel('u [m/s]');
    title(sprintf('Velocity at t = %.4f s', t_plot));

    subplot(3, 1, 3);
    plot(x_sorted, p, 'b', 'LineWidth', 2);
    grid on;
    xlabel('x [m]');
    ylabel('p [Pa]');
    title(sprintf('Pressure at t = %.4f s', t_plot))
end
