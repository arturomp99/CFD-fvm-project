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
    E_int = E - 0.5 .* rho .* u .^ 2;
    e_tot = E ./ rho; % energía total específica
    e_int = E_int ./ rho; % energía interna específica
    p = (Air.GAMMA - 1) .* (E - 0.5 .* (rhou .^ 2) ./ rho);

    [x_sorted, idx] = sort(centroids_x);
    rho = rho(idx);
    u = u(idx);
    p = p(idx);
    E = E(idx);
    E_int = E_int(idx);
    e_tot = e_tot(idx);
    e_int = e_int(idx);
    figure;

    subplot(2, 2, 1);
    plot(x_sorted, rho, 'k-', 'LineWidth', 1); hold on;
    plot(x_sorted, rho, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.6);
    grid on;
    xlabel('x [-]');
    ylabel('\rho [-]');
    title(sprintf('Density at t = %.4f [-]', t_plot));

    subplot(2, 2, 2);
    plot(x_sorted, u, 'k-', 'LineWidth', 1); hold on;
    plot(x_sorted, u, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.6);
    grid on;
    xlabel('x [-]');
    ylabel('u [-]');
    title(sprintf('Velocity at t = %.4f [-]', t_plot));

    subplot(2, 2, 3);
    plot(x_sorted, p, 'k-', 'LineWidth', 1); hold on;
    plot(x_sorted, p, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.6);
    grid on;
    xlabel('x [-]');
    ylabel('p [-]');
    title(sprintf('Pressure at t = %.4f [-]', t_plot));

    subplot(2, 2, 4);
    plot(x_sorted, e_int, 'k-', 'LineWidth', 1); hold on;
    plot(x_sorted, e_int, 'bo', 'MarkerSize', 2.5, 'LineWidth', 0.6);
    grid on;
    xlabel('x [-]');
    ylabel('e_{int} [-]');
    title(sprintf('Specific internal energy at t = %.4f [-]', t_plot));
end
